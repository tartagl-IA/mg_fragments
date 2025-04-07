"""Chemical fragment generation and filtering."""

from rdkit.Chem import BRICS, Mol, MolFromSmiles, rdMolDescriptors


def brics_from_mol(mol: Mol, min_size: int = 1) -> list[str]:
    """Generate BRICS fragments from a molecule.

    Args:
        mol (Mol): RDKit molecule object.
        min_size (int, optional): Minimum size of the BRIC fragment.
            Default is 1, meaning all fragments are returned.

    Returns:
        list[str]: List of BRICS fragments as SMILES strings.
    """
    # TODO: can be defined a minimum dimension of the BRIC with minFragmentSize
    return BRICS.BRICSDecompose(mol, minFragmentSize=min_size)


def filtered_fragments_from_mol(
    mol: Mol, min_atoms: int, max_atoms: int, flexibility: str
) -> list[Mol]:
    """Generate fragments from a molecule using BRICS decomposition and filter them based on size and flexibility.

    Args:
        mol (Mol): RDKit molecule object.
        min_atoms (int, optional): Min atoms to filter the molecule.
        max_atoms (int, optional): Max atoms to filter the molecule.
        flexibility (str, optional):
            - 'rigid' for 0 rotatable bonds
            - 'flexible' for 1 degree of freedom

    Returns:
        list[Mol]: List of filtered RDKit molecule objects.
    """
    frag_mols_list = []
    for frag in brics_from_mol(mol):
        frag_mol = MolFromSmiles(frag)
        if frag_mol and min_atoms <= frag_mol.GetNumAtoms() <= max_atoms:
            num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(frag_mol)
            if (flexibility == "rigid" and num_rotatable_bonds == 0) or (
                flexibility == "flexible" and num_rotatable_bonds == 1
            ):
                frag_mols_list.append(frag_mol)
    return frag_mols_list
