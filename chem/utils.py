"""Utilities for working with RDKit molecules."""

from io import BytesIO

from rdkit.Chem import (
    Draw,
    Mol,
    MolFromSmarts,
    MolFromSmiles,
    MolToSmarts,
    MolToSmiles,
    SDWriter,
)

from logger import get_logger

log = get_logger("Chem Utils")


def smarts_from_mol(mol: Mol) -> str:
    """Generate a SMARTS superclass string from an RDKit molecule.

    Args:
        mol (Mol): RDKit molecule object

    Returns:
        str: SMARTS string representation of the molecule
    """
    return MolToSmarts(mol)


def smiles_from_mol(mol: Mol) -> str:
    """Generate a SMILES string from an RDKit molecule.

    Args:
        mol (Mol): RDKit molecule object

    Returns:
        str: SMILES string representation of the molecule
    """
    return MolToSmiles(mol)


def mol_from_smiles(smiles: str) -> Mol:
    """Create an RDKit molecule from a SMILES string.

    Args:
        smiles (str): SMILES string representation of the molecule

    Returns:
        Mol: RDKit molecule object
    """
    return MolFromSmiles(smiles)


def mol_from_smarts(smarts: str) -> Mol:
    """Create an RDKit molecule from a SMILES string.

    Args:
        smarts (str): SMARTS string representation of the molecule

    Returns:
        Mol: RDKit molecule object
    """
    return MolFromSmarts(smarts)


def unique_mol_list(mol_list: list[Mol]) -> list[Mol]:
    """Remove duplicate molecules from a list of RDKit molecules.

    Args:
        mol_list (list[Mol]): list of RDKit molecules

    Returns:
        list[Mol]: list of unique RDKit molecules
    """
    unique_mols = {}
    for mol in mol_list:
        smi = smiles_from_mol(mol)
        if smi not in unique_mols:
            unique_mols[smi] = mol
    return list(unique_mols.values())


def mol_to_bytes(mol):
    img = Draw.MolToImage(mol, size=(300, 300))  # Generate molecule image
    buf = BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


def save_mols_to_sdf(mol_list: list[Mol], output_file: str = "fragments.sdf") -> None:
    """Save a list of RDKit molecules to an SDF file.

    Args:
        mol_list (list[Mol]): molecule list to save
        output_file (str, optional): output file name. Defaults to "fragments.sdf".
    """
    log.debug(f"Saving {len(mol_list)} molecules to SDF file: {output_file}")
    with SDWriter(output_file) as writer:
        for mol in mol_list:
            mol.SetProp("_Name", smiles_from_mol(mol))
            writer.write(mol)
    log.debug(f"Saved {len(mol_list)} molecules to SDF file: {output_file}")
