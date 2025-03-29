"""Filters for molecules."""

from rdkit.Chem import Mol, MolFromSmarts

REACTIVE_PATTERNS_MOL_LIST = [
    MolFromSmarts("C=O[C]=[C]"),  # Michael acceptor
    MolFromSmarts("[CX3]=[OX1]"),  # Aldehydes
    MolFromSmarts("[CX3](=O)[Cl,Br,I]"),  # Acyl halides
    MolFromSmarts("[CX3](=O)N=C=O"),  # Isocyanates
    MolFromSmarts("[SX4](=O)(=O)[Cl,Br,I]"),  # Sulfonyl halides
    MolFromSmarts("[C;R][O;R][C;R]"),  # Epoxides
    MolFromSmarts("[C;R][N;R][C;R]"),  # Aziridines
    MolFromSmarts("[C;R][S;R][C;R]"),  # Thiiranes
]


def mol_is_reactive(mol: Mol) -> bool:
    """Check if a molecule contains electrophilic groups using precompiled SMARTS.

    Args:
        mol (Mol): Molecule to check.

    Returns:
        bool: True if the molecule contains reactive groups, False otherwise.
    """
    return any(mol.HasSubstructMatch(pattern) for pattern in REACTIVE_PATTERNS_MOL_LIST)
