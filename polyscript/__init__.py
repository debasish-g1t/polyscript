"""PolyScript: A chemical notation format and toolkit for polymer sequences.

PolyScript encodes polymer synthesis information as a structured string
containing monomers (SMILES), a reaction (SMARTS), and the resulting
polymer (SMILES/PSMILES), e.g.::

    ${monomer1}+{monomer2}=>{reaction_smarts}=>{polymer}$

The package provides utilities for parsing, validating, converting, and
analysing PolyScript sequences, including MPI-parallel batch validation
of large datasets.
"""
