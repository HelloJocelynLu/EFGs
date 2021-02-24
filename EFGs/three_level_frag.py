import itertools
import re
import argparse
from collections import deque
from itertools import combinations

import rdkit
import rdkit.Chem as Chem

from .ifg import identify_functional_groups

lg = rdkit.RDLogger.logger() 
lg.setLevel(rdkit.RDLogger.CRITICAL)

patt = r'[C,H][0-9]{2}[0,-1,1]'
IsH010 = lambda a: (a.GetSymbol()=='H' and a.GetDegree()==1 and a.GetFormalCharge()==0)
IsAromatic = lambda a: a.GetIsAromatic()

class WordNotFoundError(Exception):
    pass

def standize(smiles, RemoveMap=True, canonical=True, isomericSmiles=True, Order=False, asMol=False):
    '''    standize(smiles, RemoveMap=True, canonical=True, isomericSmiles=True, Order=False, asMol=False)
    Parameters
    ----------
        smiles: molecular SMILES representations or rdkit.Chem.rdchem.Mol (if asMol=True)
            Molecule
        RemoveMap: bool (default: True)
            Remove atomic mapped number
        canonical: bool (default: True)
            Canonicalize
        isomericSmiles: bool (default=True)
            Take chirality into considerations.
        Order: bool (default: False)
            Return original atom indices based on output SMILES order 
        asMol: bool (default: False)
            Input molecule is an rdkit.Chem.rdchem.Mol object
    Returns
    -------
        standize_SMILES: str
            Standized SMILES based on input rules
        order: tuple (optional)
            Return only when Order=True'''

    if not smiles:
        return ''
    if asMol:
        mol = smiles.__copy__()
    else:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if not mol:
        return ''
    if not RemoveMap:
        return Chem.MolToSmiles(mol, canonical=canonical)
    try:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        raw_smiles = Chem.MolToSmiles(mol, canonical=canonical)
        init_order = list(mol.GetPropsAsDict(True,True)['_smilesAtomOutputOrder'])
        # RDKit bug: sometimes only one canonicalization is not enough, e.g.: 'Nc1cc2c3ccc(n1)c2c3' when asMol
        mol2 = Chem.MolFromSmiles(raw_smiles)
        sanitized = Chem.MolToSmiles(mol2, canonical=canonical, isomericSmiles=isomericSmiles)
        next_order = list(mol2.GetPropsAsDict(True,True)['_smilesAtomOutputOrder'])
        if not Order:
            return sanitized
        match = tuple(init_order[i] for i in next_order)
        return sanitized, match
    except:
        return ''

def sp3merge(atom1, atom2):
    '''Given two atoms, to see if they can merge into a simple carbon chain.'''
    if atom1.GetSymbol() == atom2.GetSymbol() == 'C':
        if atom1.GetIsAromatic() == atom2.GetIsAromatic() == False and atom1.GetDegree()<3 and atom2.GetDegree()<3:
            return Chem.BondType.SINGLE
    return False

def aromaticmerge(atom1, atom2, bond):
    '''Given two atoms and bond between them, to see if they can merge into an aromatic ring.'''
    # I would not expect structures like `c1ccc2c(c1)-c1ccccc1-2`
    if atom1.GetIsAromatic() and atom2.GetIsAromatic() and bond.GetIsAromatic():
        return Chem.BondType.AROMATIC
    # Also include some dipoles
    if atom1.GetIsAromatic() or atom2.GetIsAromatic():
        if bond.GetBondTypeAsDouble()==2.0 or ((atom1.GetFormalCharge()!=0 and atom1.GetIsAromatic()) or \
        (atom2.GetFormalCharge()!=0 and atom2.GetIsAromatic())):
            return Chem.BondType.DOUBLE
    return False

def AtomListToSubMol(mol, amap, bmap=(), includeConformer=False):
    """AtomListToSubMol(mol, amap, bmap=(), includeConformer=False)
    Parameters
    ----------
        mol: rdkit.Chem.rdchem.Mol
            Molecule
        amap: array-like
            List of atom indices (zero-based)
        bmap: array-like
            List of bond indices (zero-based. Defult=range(mol.GetNumBonds()))
        includeConformer: bool (default=False)
            Toogle to include atoms coordinates in submolecule.

    Returns
    -------
        submol: rdkit.Chem.rdchem.RWMol
            Submol determined by specified atom list
    """
    if not bmap: bmap=range(mol.GetNumBonds())
    if not isinstance(amap, list):
        amap = list(amap)
    if not isinstance(bmap, list):
        bmap = list(bmap)
    submol = Chem.RWMol(Chem.MolFromSmiles(''))
    for aix in amap:
        submol.AddAtom(mol.GetAtomWithIdx(aix))
    for i, j in combinations(amap, 2):
        bond = mol.GetBondBetweenAtoms(i, j)
        if bond and bond.GetIdx() in bmap:
            submol.AddBond(amap.index(i),
                           amap.index(j),
                           bond.GetBondType())
    if includeConformer:
        for conf in mol.GetConformers():
            new_conf = Chem.Conformer(len(amap))
            for i in range(len(amap)):
                new_conf.SetAtomPosition(i, conf.GetAtomPosition(amap[i]))
                new_conf.SetId(conf.GetId())
                new_conf.Set3D(conf.Is3D())
            submol.AddConformer(new_conf)
    return submol

def merge(iterable, style=set):
    LL = set(itertools.chain.from_iterable(iterable))
    for ele in LL:
        components = [x for x in iterable if ele in x]
        for i in components:
            iterable.remove(i)
        iterable += [style(set(itertools.chain.from_iterable(components)))]

def extractAromatic(mol):
    '''Given a mol object, return the 'aromatic' part in that compound.
    Return: aro: a list of smiles of aromatic structures
    idx: a list of atom indexes corresponding to aromatic structures'''
    def san_check(mole, atom_maps):
        left = set(range(mole.GetNumAtoms())).difference(atom_maps)
        if any([mole.GetAtomWithIdx(i).GetIsAromatic() for i in left]):
            atom_maps.extend([i for i in left if mole.GetAtomWithIdx(i).GetIsAromatic()])
            raise WordNotFoundError
    def NotSingleBond(mole, atom1, atom2):
        bond = mole.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx())
        if not bond:
            return False
        if bond.GetBondType() == Chem.BondType.SINGLE:
            return False
        return True
    m = mol.__copy__()
    for atom in m.GetAtoms():
        atom.SetAtomMapNum(0)
    aro, new_mol_index = [], []
    amap, bmap=set(), set()
    Chem.Kekulize(m)
    rinfo = m.GetRingInfo()
    n_atom = rinfo.AtomRings()
    for bond in m.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        order = aromaticmerge(atom1, atom2, bond)
        if order:
            amap.update([atom1.GetIdx(), atom2.GetIdx()])
            bmap.add(bond.GetIdx())
    amap = list(amap)
    bmap = list(bmap)
    rwmol = AtomListToSubMol(m, amap, bmap)
    try:
        san_check(m, amap)
        Chem.SanitizeMol(rwmol)
    except:
        rwmol = AtomListToSubMol(m, amap)
    groups = []
    for mm in Chem.GetMolFrags(rwmol):
        cur_mol = AtomListToSubMol(rwmol, mm)
        std_smiles = standize(Chem.MolToSmiles(cur_mol))
        cur_idx = [amap[i] for i in mm]
        if any([atom.GetIsAromatic() for atom in Chem.MolFromSmiles(std_smiles).GetAtoms()]) and \
            [IsAromatic(mol.GetAtomWithIdx(i)) for i in cur_idx] == [IsAromatic(rwmol.GetAtomWithIdx(i)) for i in mm]:
            updated = set(cur_idx)
        else:
            updated = set([x for group in n_atom for x in group if len(set(cur_idx)&set(group))>2])
        queue = deque(updated)
        while len(queue) > 0:
            extra_atom = queue.popleft()
            EA = m.GetAtomWithIdx(extra_atom)
            if set([aa.GetIdx() for aa in EA.GetNeighbors() if NotSingleBond(m, EA, aa)])<=updated:
                continue
            arm = set([aa.GetIdx() for aa in EA.GetNeighbors() if NotSingleBond(m, EA, aa)]).difference(updated)
            queue.extend(arm)
            updated.update(arm)
        groups.append(updated)
    merge(groups)
    Chem.SanitizeMol(m)
    Chem.Kekulize(m)
    for group in groups:
        group = list(group)
        cur_mol = AtomListToSubMol(m, group)
       # Chem.SanitizeMol(cur_mol)
        aro_fg, order = standize(cur_mol, Order=True, asMol=True)
        aro.append(standize(aro_fg))
        new_mol_index.append(tuple(group[i] for i in order))
    return aro, new_mol_index

def aro_ifg(mol, idx2map=(), mergeAtom=True, isomericSmiles=True):
    aro, aro_idx = extractAromatic(mol)
    if idx2map:
        masked = [idx2map[x] for i in aro_idx for x in i]
        aro_idx = [tuple(map(lambda t:idx2map[t],x)) for x in aro_idx]
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(idx2map[atom.GetIdx()])
        UseMap = True
    else:
        masked = [x for i in aro_idx for x in i]
        UseMap = False
    fgs = identify_functional_groups(mol, MapNum=UseMap, masked=set(masked), mergeAtom=mergeAtom, isomericSmiles=isomericSmiles)
    return aro, aro_idx, fgs

def breakBond(m, MapNum=False, returnidx=False):
    '''Given a mol, this function would break all the single bond in the molecule,
    and then return (fragments smiles, Carbons/hydrogens Types)'''
    mol = m.__copy__()
    idx2map = {}
    CHs, CHs_idx, frags, frag_idx = [], [], [], []
    for atom in mol.GetAtoms():
        if MapNum:
            idx2map[atom.GetIdx()]=atom.GetAtomMapNum()
        atom.SetAtomMapNum(atom.GetIdx())
    no_link = set()
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() == 1.0:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if sp3merge(atom1, atom2):
                continue
            no_link.add(bond.GetIdx())
    links = set(range(mol.GetNumBonds())).difference(no_link)
    if not links: links=[mol.GetNumBonds()]
    # To avoid empty sp3_link ==> all bonds, set an impossible value
    rwmol = AtomListToSubMol(mol, range(mol.GetNumAtoms()), links)
    for mm in Chem.GetMolFrags(rwmol, asMols=True, sanitizeFrags=False):
        if mm.GetNumAtoms() == 1 and mm.GetAtomWithIdx(0).GetSymbol() in ['C', 'H']:
            atom = mol.GetAtomWithIdx(mm.GetAtomWithIdx(0).GetAtomMapNum())
            CHs.append(atom.GetSymbol()+str(int(atom.GetIsAromatic()))+str(atom.GetDegree())
        +str(atom.GetFormalCharge()))
            CHs_idx.append((atom.GetAtomMapNum(),))
        else:
            mm_renum, order = standize(mm, asMol=True, Order=True)
            frag_idx.append(tuple([mm.GetAtomWithIdx(i).GetAtomMapNum() for i in order]))
            frags.append(mm_renum)
    if not returnidx:
        return frags, CHs
    if idx2map:
        frag_idx = [tuple(map(lambda x: idx2map[x],group)) for group in frag_idx]
        CHs_idx = [tuple(map(lambda x: idx2map[x],group)) for group in CHs_idx]
    return frags, CHs, frag_idx, CHs_idx

def mol2frag(raw_mol, TreatHs='ignore', returnidx=False, toEnd=False, vocabulary=(), extra_included=False, isomericSmiles=True, UnknownIdentity=False, extra_backup={}):
    '''
    raw_mol: rdkit mol object to be decompose
    TreatHs: (optional) The way to treat Hs. Default: 'ignore' (Other options: 'separate': treat Hs separately; 'include': merged to neighboring EFGs) 
    returnidx: (optional) Whether cooresponding atom index of EFGs shouls be returned. (Default=False)
    toEnd: Whether to decompose to the end. (Default=False, will only do 1-step decomposition)
    vocabulary: (optional) A list of smiles which contains EFGs. This argument would be ignore if 
    toEnd=False. If toEnd is set to True, this argument is required. (Default=None)
    extra_included: (optional) If fragments outside of vocabulary should be parsed. (Default=False, will throw
    an error if a fragment cannot be found in vocabulary). When it is set to True, additional fragments
    would be simply classified based on their aromaticity.
    isomericSmiles: (optional) include information about stereochemistry in the SMILES.
    extra_backup: (optional) If an empty dictionary is provided, additional fragments' smiles would be added. 
    return:
    Functional groups' smiles (or 'Aromatic'/'Non-aromatic') and C/Hs
    (or)
    Functional groups' smiles (or 'Aromatic'/'Non-aromatic'), C/Hs, atom indices of funtional groups
    and CHs (if returnidx=True)
    '''
    mol = raw_mol.__copy__()
    CHs = []
    CHs_idx = []
    idx2map = {}
    # Primary decomposition
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
        idx2map[atom.GetIdx()]=atom.GetAtomMapNum()
    if TreatHs != 'ignore':
        if TreatHs == 'include':
            h_bond = []
            for bond in mol.GetBonds():
                atom1 = bond.GetBeginAtom()
                atom2 = bond.GetEndAtom()
                if any([atom.GetSymbol()=='H' for atom in (atom1, atom2)]):
                    h_bond.append((atom1.GetIdx(), atom2.GetIdx()))
        mol = Chem.RemoveHs(mol)
        non_H = [atom.GetAtomMapNum() for atom in mol.GetAtoms()]
    more_frag, frag_idx, fgs = aro_ifg(mol, idx2map=idx2map, isomericSmiles=isomericSmiles)
    # Further decomposition
    if toEnd:
        if len(vocabulary)==0: raise WordNotFoundError('No volcabulary but toEnd is set to True.')
        if extra_included:
            if not 'Aromatic' in vocabulary:
                vocabulary.append('Aromatic')
            if not 'Non-aromatic' in vocabulary:
                vocabulary.append('Non-aromatic')
        else:
            if (('Aromatic' in vocabulary) or ('Non-aromatic' in vocabulary)) and not extra_included:
                raise WordNotFoundError('Please either remove wild cards ("Aromatic"/"Non-aromatic") in vocabulary or set "extra_included=True"')
        # If primary decomposition does not work, simply do 'shallow' merge
        if not set([fg.atoms for fg in fgs])<=set(vocabulary):
            extra = [fg for fg in fgs if fg.atoms not in vocabulary]
            White_list = []
            for fg in extra:
                fg_smiles = fg.atoms
                fg_idx = fg.atomIds
                submol = smiles2mol(fg_smiles, fg_idx)
                temp_aro, temp_aro_idx, temp_frags = aro_ifg(Chem.MolFromSmiles(fg_smiles), idx2map={atom.GetIdx():atom.GetAtomMapNum() for atom in submol.GetAtoms()}, mergeAtom=False, isomericSmiles=isomericSmiles)
                if any([Chem.MolFromSmiles(x)==None for x in [temp_fg.atoms for temp_fg in temp_frags]+temp_aro]):
                    White_list.append(fg)
                    continue
                else:
                    fgs.remove(fg)
                    more_frag.extend(temp_aro)
                    frag_idx.extend(temp_aro_idx)
                    fgs.extend(temp_frags)
                    for temp_fg in temp_frags:
                        temp_mol = smiles2mol(temp_fg.atoms, temp_fg.atomIds)
                        if any([atom.GetIsAromatic() for atom in temp_mol.GetAtoms()]):
                            a,b,c = aro_ifg(Chem.MolFromSmiles(temp_fg.atoms), idx2map={atom.GetIdx():atom.GetAtomMapNum() for atom in temp_mol.GetAtoms()}, mergeAtom=False, isomericSmiles=isomericSmiles)
                            if any([Chem.MolFromSmiles(x)==None for x in a+[m.atoms for m in c]]):
                                White_list.append(temp_fg)
                                continue
                            else:
                                more_frag.extend(a)
                                frag_idx.extend(b)
                                fgs.remove(temp_fg)
                                fgs.extend(c)
            if not set([fg.atoms for fg in fgs])<=set(vocabulary):
                extra = [fg for fg in fgs if fg.atoms not in vocabulary]
                for fg in extra:
                    if fg in White_list: continue
                    submol = smiles2mol(fg.atoms, fg.atomIds)
                    try:
                        a, b, c, d = breakBond(submol, MapNum=True, returnidx=True)
                        fgs.remove(fg)
                        more_frag.extend(a)
                        frag_idx.extend(c)
                    # Fragments should be valid
                    except argparse.ArgumentError:
                        pass
    fgs_atoms = [num for ids in fgs for num in ids.atomIds] + [num2 for tuples in frag_idx for num2 in tuples]
    fgs_tuple = [ids.atomIds for ids in fgs] + [tuples for tuples in frag_idx]
    n_atom = mol.GetNumAtoms()
    # to see if we have any atoms left behind
    rest_atom = {i for i in range(n_atom)}.difference(fgs_atoms)
    sp3_link = set()
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if set((atom1.GetIdx(), atom2.GetIdx()))<=rest_atom and sp3merge(atom1, atom2):
            sp3_link.add(bond.GetIdx())
    if not sp3_link: sp3_link=[mol.GetNumBonds()]
    # To avoid empty sp3_link ==> all bonds, set an impossible value
    rwmol = AtomListToSubMol(mol, rest_atom, sp3_link)
    for mm in Chem.GetMolFrags(rwmol, asMols=True, sanitizeFrags=False):
        if mm.GetNumAtoms() == 1 and mm.GetAtomWithIdx(0).GetSymbol() in ['C', 'H']:
            atom = mol.GetAtomWithIdx(mm.GetAtomWithIdx(0).GetAtomMapNum())
            if TreatHs == 'include' and IsH010(atom): continue
            CHs.append(atom.GetSymbol()+str(int(atom.GetIsAromatic()))+str(atom.GetDegree())
        +str(atom.GetFormalCharge()))
            CHs_idx.append((atom.GetAtomMapNum(),))
        else:
            mm_renum, order = standize(mm, asMol=True, Order=True)
            CHs_idx.append(tuple([mm.GetAtomWithIdx(i).GetAtomMapNum() for i in order]))
            CHs.append(mm_renum)
    nonCHs=[fg.atoms for fg in fgs]+more_frag
    if not set(nonCHs+CHs)<=set(vocabulary) and len(vocabulary)>0:
        if not extra_included:
            raise WordNotFoundError('{} not found.'.format(set(nonCHs+CHs).difference(vocabulary)))
        else:
            extras = set(nonCHs+CHs).difference(vocabulary)
            for extra in extras:
                append_context = ''
                if UnknownIdentity:
                    append_context = '.'+extra
                if re.match(patt, extra):
                    CHs = ['Non-aromatic'+append_context if x==extra else x for x in CHs]
                    try:
                        extra_backup['Non-aromatic'].add(extra)
                    except KeyError:
                        extra_backup['Non-aromatic']={extra}
                    continue
                sub_mol = Chem.MolFromSmiles(extra)
                if any([atom.GetIsAromatic() for atom in sub_mol.GetAtoms()]):
                    nonCHs = ['Aromatic'+append_context if x==extra else x for x in nonCHs]
                    CHs = ['Aromatic'+append_context if x==extra else x for x in CHs]
                    try:
                        extra_backup['Aromatic'].add(extra)
                    except KeyError:
                        extra_backup['Aromatic']={extra}
                else:
                    nonCHs = ['Non-aromatic'+append_context if x==extra else x for x in nonCHs]
                    CHs = ['Non-aromatic'+append_context if x==extra else x for x in CHs]
                    try:
                        extra_backup['Non-aromatic'].add(extra)
                    except KeyError:
                        extra_backup['Non-aromatic']={extra}
    if TreatHs !='ignore':
        H_s = set(range(len(idx2map))).difference(non_H)
        CHs_idx = [tuple(map(lambda t:non_H[t], i)) for i in CHs_idx]
        fgs_tuple = [tuple(map(lambda t:non_H[t], i)) for i in fgs_tuple]
        if TreatHs == 'include':
            for pos in range(len(CHs_idx)):
                cur_group = set(CHs_idx[pos])
                for bond in list(h_bond):
                    if cur_group&set(bond):
                        cur_group.update(bond)
                        h_bond.remove(bond)
                CHs_idx[pos] = tuple(cur_group)
            for pos in range(len(fgs_tuple)):
                cur_group = set(fgs_tuple[pos])
                for bond in list(h_bond):
                    if cur_group&set(bond):
                        cur_group.update(bond)
                        h_bond.remove(bond)
                fgs_tuple[pos] = tuple(cur_group)
            H_s = set(itertools.chain.from_iterable(h_bond))
            # H_s is not empty only for H2
        for hid in H_s:
            CHs.append('H010')
            CHs_idx.append((hid,))
    if not returnidx:
        return nonCHs, CHs
    return nonCHs, CHs, fgs_tuple, CHs_idx

def counter(x, dic, increase = 1):
    '''A counter to record the occurance of x in dic'''
    if x in dic:
        dic[x] += increase
    else:
        dic[x] = increase

def smiles2mol(smiles, maplist=None):
    '''Given a smiles, return a mol'''
    mol = Chem.MolFromSmiles(smiles)
    if not mol: mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol.UpdatePropertyCache()
    if maplist:
        for i,atom in enumerate(mol.GetAtoms()):
            atom.SetAtomMapNum(maplist[i])
    return mol

def _cleavage1(dictionary, White_list, alpha = 0.7, isomericSmiles=True):
    patt = r'[C,H][0-9]{2}[0,-1,1,4]{0,1}'
    qua = sorted([dictionary[x] for x in dictionary])[::-1][int(alpha*len(dictionary))-1]
    rare = [x for x in dictionary if dictionary[x]<=qua]
    for smi in rare:
        if re.match(patt, smi) or (smi in White_list): continue
        num = dictionary[smi]
        del dictionary[smi]
        mol = Chem.MolFromSmiles(smi)
        temp_aro, temp_aro_idx, temp_frags = aro_ifg(mol, mergeAtom=False, isomericSmiles=isomericSmiles)
        # The whole molecule is an aromatic species
        if not temp_frags:
            White_list.append(smi)
            counter(smi, dictionary, increase=num)
            continue
        frags = [temp_fg.atoms for temp_fg in temp_frags]+temp_aro
        if any([Chem.MolFromSmiles(x)==None for x in frags]):
            White_list.append(smi)
            counter(smi, dictionary, increase=num)
            continue
        extracted = [ids.atomIds for ids in temp_frags]+temp_aro_idx
        rest_atom = set(range(mol.GetNumAtoms())).difference([num for ids in extracted for num in ids])
        rwmol = Chem.RWMol(Chem.MolFromSmiles(''))
        CHs = []
        passed = {}
        for bond in mol.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if set((atom1.GetIdx(), atom2.GetIdx()))<rest_atom:
                if atom1.GetIdx() not in passed:
                    passed[atom1.GetIdx()] = rwmol.AddAtom(atom1)
                if atom2.GetIdx() not in passed:
                    passed[atom2.GetIdx()] = rwmol.AddAtom(atom2)
                order = sp3merge(atom1, atom2)
                if order:
                    rwmol.AddBond(passed[atom1.GetIdx()], passed[atom2.GetIdx()], order=bond.GetBondType())
        passed_ = {value:key for key, value in passed.items()}
        new_mol_index = [tuple(map(lambda t:passed_[t], i)) for i in Chem.GetMolFrags(rwmol)]
        for atom in rwmol.GetAtoms():
            atom.SetAtomMapNum(0)
        for mid, mm in zip(new_mol_index, Chem.GetMolFrags(rwmol, asMols=True, sanitizeFrags=False)):
            if mm.GetNumAtoms() == 1 and mm.GetAtomWithIdx(0).GetSymbol() in ['C', 'H']:
                atom = mol.GetAtomWithIdx(mid[0])
                CHs.append(atom.GetSymbol()+str(int(atom.GetIsAromatic()))+str(atom.GetDegree())
            +str(atom.GetFormalCharge()))
            else:
                CHs.append(Chem.MolToSmiles(mm))
        for fg in frags:
            counter(fg, dictionary, increase=num)
        for ch in CHs:
            counter(ch, dictionary, increase=num)
            
def _cleavage2(dictionary, White_list, alpha = 0.7, isomericSmiles=True):
    patt = r'[C,H][0-9]{2}[0,-1,1,4]{0,1}'
    qua = sorted([dictionary[x] for x in dictionary])[::-1][int(alpha*len(dictionary))-1]
    rare = [x for x in dictionary if dictionary[x]<=qua]
    for smi in rare:
        if re.match(patt, smi) or (smi in White_list): continue
        num = dictionary[smi]
        del dictionary[smi]
        mol = Chem.MolFromSmiles(smi)
        temp_aro, temp_aro_idx, temp_frags = aro_ifg(mol, mergeAtom=False, isomericSmiles=isomericSmiles)
        # The whole molecule is an aromatic species
        if not temp_frags:
            White_list.append(smi)
            counter(smi, dictionary, increase=num)
            continue
        try:
            a, b, c, d = breakBond(mol, returnidx=True)
            # Fragments should be valid
            White_list.append(smi)
            counter(smi, dictionary, increase=num)
            continue
        except argparse.ArgumentError:
            pass
        for fg in a:
            counter(fg, dictionary, increase=num)
        for ch in b:
            counter(ch, dictionary, increase=num)

def cleavage(dictionary, alpha = 0.7, isomericSmiles=True):
    old_size = len(dictionary)
    white_list=[]
    _cleavage1(dictionary, white_list, alpha = alpha, isomericSmiles=isomericSmiles)
    new_size = len(dictionary)
    while old_size-new_size>0:
        old_size = len(dictionary)
        _cleavage1(dictionary, white_list, alpha = alpha, isomericSmiles=isomericSmiles)
        white_list = list(set(white_list))
        new_size = len(dictionary)
    _cleavage2(dictionary, white_list, alpha = alpha, isomericSmiles=isomericSmiles)
    new_size = len(dictionary)
    while old_size-new_size>0:
        old_size = len(dictionary)
        _cleavage2(dictionary, white_list, alpha = alpha, isomericSmiles=isomericSmiles)
        white_list = list(set(white_list))
        new_size = len(dictionary)
    # To remove rare EFGs
    qua = sorted([dictionary[x] for x in dictionary])[::-1][int(alpha*len(dictionary))-1]
    rare = [x for x in dictionary if dictionary[x]<=qua]
    for x in rare:
        del dictionary[x]
