EFGs (Extended functional groups)
====

.. image:: https://img.shields.io/pypi/v/EFGs.svg
    :target: https://pypi.python.org/pypi/EFGs
    :alt: Latest PyPI version

.. image:: https://travis-ci.org/borntyping/cookiecutter-pypackage-minimal.png
   :target: https://travis-ci.org/borntyping/cookiecutter-pypackage-minimal
   :alt: Latest Travis CI build status

Extended Functional Groups

Extended functional group is a generalized version of traditional functional group and it also contains chemical groups that formed by only carbon atoms. It is inspired by `Peter Ertl`_'s work: 

Ertl, P. An algorithm to identify functional groups in organic molecules. *J Cheminform* **9**, 36 (2017)

.. _Peter Ertl: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0225-z 

Built based on that, we also induced the idea that a moelcule should be fully covered by 'Functional Groups'.

The philosophy of EFG (Extended functional group) is to do fragmentation on molecules so that all fragments of the molecule are chemical valid. To do that, we:

1. **Identify aromatic structures.** If two atoms shared the same aromatic ring system, they would be merged.
2. **Identify special substructures**:
    * Mark all heteroatoms in a molecule
    * Mark ‘special’ carbon atoms (carbon atoms with double/triple bonds, acetal carbons and three-membered heterocycles.)
    * Merge all connected marked atoms to a single functional group
3. **Identify simple carbon chains**: sp3 carbons connected by two or more hydrogens
4. **Other single atoms** The number of single atoms can be significantly reduced by defining subclasses and merging some of them together. All atoms are classified by their aromaticity, degree and formal charge and recorded as element symbol followed by three number corresponding to above properties. For example, Hydrogen ($𝐻_2$) would be H010, methyl group would be C010.

.. image:: image.png

In order to alleviate the imbalance distribution of different EFGs, we proposed an iterative way to selectively decompose large functional groups:

1. Set a cut-off value α (0<α<1)

2. Collect sparse functional groups whose rankings are behind top α in frequency distribution

3. Further decompose collected functional groups:

    * a. Neighboring small functional groups which would be merged before would not be merged anymore unless they have shared atom(s).
    * b. (If i. is not applicable) Cut all single bonds
4. Repeat previous steps until the number of functional groups does not change.

For most molecular datasets, this method is able to describe > 99% molecules with < 1% number of EFGs. 

Usage
-----

See *Tutorial.ipynb* in Examples/ folder for detailed examples.

*mol2frag* is the core function to do the fragmentation.

Installation
------------

.. code:: bash

   $ python setup.py install

Requirements
^^^^^^^^^^^^

rdkit >= 2019.03

Licence
-------

Authors
-------

`EFGs` was written by `Jocelyn Lu <jl8570@nyu.edu>`_.
