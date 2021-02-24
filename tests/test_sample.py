# Sample Test passing with nose and pytest
import pytest
from rdkit import Chem
from EFGs import mol2frag

def test_standize():
    from EFGs import standize
    with open("tests/custom_tests.txt") as rf:
        for line in rf:
            standize(line.strip())
    assert True

def test_pass():
    with open("tests/custom_tests.txt") as rf:
        for line in rf:
            smiles = line.strip()
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                raise ValueError("invalid input:{}".format(smiles))
            nonCH, CH, fgidx, chidx = mol2frag(mol, returnidx=True)
    assert True, "Test failed!"

