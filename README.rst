EFGs
====

.. image:: https://img.shields.io/pypi/v/EFGs.svg
    :target: https://pypi.python.org/pypi/EFGs
    :alt: Latest PyPI version

.. image:: https://travis-ci.org/borntyping/cookiecutter-pypackage-minimal.png
   :target: https://travis-ci.org/borntyping/cookiecutter-pypackage-minimal
   :alt: Latest Travis CI build status

Extended Functional Groups

Version 0.1.0
Original Version


Version 0.2.0
1. Identify functional group: atomID and typeID are arranged based on index order in atom and type.
2. mol2frag, modified the code (remove the calling of GetSubstructMatches)
3. mol2frag, fix some substructure decomposition with 'white list'

Version 0.3.0
1. mol2frag: And new argument: isomericSmiles=True
2. mol2frag: Fix the bug that got wrong index (indexes changed even if canonical is set to False) when deal with aromatic complex molecules 
3. cleavage, fix the iterative algorithm

Version 0.3.4
1. standize: accept mol as inputs. Be able to return order if Order=True
2. Fix bugs when dealing with some molecules.
e.g.:
N#Cc1cccc2[nH][nH]c3nc(c1)-c23
O=C=c1cc2ccc3c(c2[nH]1)C=CN=3

Version 0.4.0
1. extractAromatic: Improved with merge algotirhm. Fix bugs when dealing with some molecules.
e.g.:
N#Cc1ccc2c(c1)NC(Cl)=c1c(Cl)ccnc1=N2

Usage
-----

Installation
------------

Requirements
^^^^^^^^^^^^

Compatibility
-------------

Licence
-------

Authors
-------

`EFGs` was written by `Jocelyn Lu <jl8570@nyu.edu>`_.
