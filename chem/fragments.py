"""Chemical fragment generation and filtering."""

from rdkit.Chem import BRICS, Mol, MolFromSmiles, rdMolDescriptors


def filtered_fragments_from_mol(
    mol: Mol, min_atoms: int = 5, max_atoms: int = 100, flexibility: str = "rigid"
) -> list[Mol]:
    """Generate fragments from a molecule using BRICS decomposition and filter them based on size and flexibility.

    Args:
        mol (Mol): RDKit molecule object.
        min_atoms (int, optional): Min atoms to filter the molecule. Defaults to 5.
        max_atoms (int, optional): Max atoms to filter the molecule. Defaults to 100.
        flexibility (str, optional):
            - 'rigid' for 0 rotatable bonds
            - 'flexible' for 1 degree of freedom. Defaults to "rigid"

    Returns:
        list[Mol]: List of filtered RDKit molecule objects.
    """
    frag_mols_list = []
    for frag in BRICS.BRICSDecompose(mol):
        frag_mol = MolFromSmiles(frag)
        if frag_mol and min_atoms <= frag_mol.GetNumAtoms() <= max_atoms:
            num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(frag_mol)
            if (flexibility == "rigid" and num_rotatable_bonds == 0) or (
                flexibility == "flexible" and num_rotatable_bonds == 1
            ):
                frag_mols_list.append(frag_mol)
    return frag_mols_list
