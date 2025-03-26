from rdkit.Chem import Mol, MolFromSmarts, MolFromSmiles


ELECTROPHILIC_PATTERNS = [
    MolFromSmarts("C=O[C]=[C]"),                # Michael acceptor
    MolFromSmarts("[CX3]=[OX1]"),               # Aldehydes
    MolFromSmarts("[CX3](=O)[Cl,Br,I]"),        # Acyl halides
    MolFromSmarts("[CX3](=O)N=C=O"),            # Isocyanates
    MolFromSmarts("[SX4](=O)(=O)[Cl,Br,I]"),    # Sulfonyl halides
    MolFromSmarts("[C;R][O;R][C;R]"),           # Epoxides
    MolFromSmarts("[C;R][N;R][C;R]"),           # Aziridines
    MolFromSmarts("[C;R][S;R][C;R]"),           # Thiiranes
]

def mol_from_smiles(smiles: str) -> Mol:
    """
    Create an RDKit molecule from a SMILES string.
    """
    return MolFromSmiles(smiles)


def mol_is_reactive(mol: Mol) -> bool:
    """
    Check if a molecule contains electrophilic groups using precompiled SMARTS.
    """
    return any(mol.HasSubstructMatch(pattern) for pattern in ELECTROPHILIC_PATTERNS)
