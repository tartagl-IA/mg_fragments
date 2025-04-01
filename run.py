"""MG Fragments runner.

This script processes a list of target IDs, retrieves their associated molecules
from a database, filters the fragments of these molecules, and saves the unique fragments
to SDF files.
"""

import os
import sys

import settings
from chem.filters import mol_reactive
from chem.fragments import Mol, filtered_fragments_from_mol
from chem.utils import (
    mol_from_smiles,
    save_mols_to_sdf,
    smarts_from_mol,
    unique_mol_list,
)
from db_mg_fragments import close_db_connection as db_mgf_close_db_connection
from db_mg_fragments import get_db_connection as db_mgf_get_db_connection
from db_mg_fragments.handlers import mols as db_mgf_mols_handlers
from logger import get_logger

log = get_logger(sys.modules[__name__].__name__)


def filtered_fragments_from_smiles(smiles: str) -> None | dict[str, list[Mol]]:
    """Get filtered fragments from a SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        None | dict[str, list[Mol]]: Dictionary of SMARTS and associated Molecules
        or None if the molecule is reactive.
    """
    mol = mol_from_smiles(smiles)
    if mol_reactive(mol):
        log.debug(f"Mol is reactive: {smiles}")
        return None
    fragment_mol_list = filtered_fragments_from_mol(
        mol,
        settings.FRAGMENTS_MIN_ATOMS,
        settings.FRAGMENTS_MAX_ATOMS,
        settings.FRAGMENTS_FLEXIBILITY,
    )

    smarts_fragments_dict: dict[str, list[Mol]] = {}
    for fragment_mol in fragment_mol_list:
        smarts = smarts_from_mol(fragment_mol)
        if smarts not in smarts_fragments_dict:
            smarts_fragments_dict[smarts] = []
        smarts_fragments_dict[smarts].append(fragment_mol)
    return smarts_fragments_dict


def get_top_results(
    smarts_fragments_dict: dict[str, list[Mol]], top_n: int
) -> dict[str, list[Mol]]:
    """Get the top N results from a dictionary of SMILES fragments.

    Args:
        smarts_fragments_dict (dict[str, list[Mol]]): dictionary of SMARTS and associated Molecules
        top_n (int, optional): Retrieve first N results.

    Returns:
        dict[str, list[Mol]]: dictionary of top N SMARTS and associated Molecules
    """
    smarts_fragments_dict_sorted = dict(
        sorted(
            smarts_fragments_dict.items(), key=lambda item: len(item[1]), reverse=True
        )
    )
    return dict(list(smarts_fragments_dict_sorted.items())[:top_n])


def main(target_id_list: list[str]) -> None:
    """Main function to process target IDs and save fragments to SDF files."""
    db_mgf_connection = db_mgf_get_db_connection()
    for target_id in target_id_list:
        log.info(f"Processing target: {target_id}")
        target_smarts_fragments_dict = dict()
        for mol_db in db_mgf_mols_handlers.get_by_target(db_mgf_connection, target_id):
            smarts_fragments_dict = filtered_fragments_from_smiles(
                mol_db["canonical_smiles"]
            )
            if smarts_fragments_dict is None:
                continue
            target_smarts_fragments_dict.update(smarts_fragments_dict)
        smarts_fragments_dict_top = get_top_results(
            target_smarts_fragments_dict, settings.FRAGMENTS_TOP_RES
        )
        fragments_mol_top_all = []
        for mol_list in smarts_fragments_dict_top.values():
            fragments_mol_top_all.extend(mol_list)
        unique_fragment_mols = unique_mol_list(fragments_mol_top_all)
        if not os.path.exists(settings.FRAGMENTS_OUTPUT_DIR):
            os.makedirs(settings.FRAGMENTS_OUTPUT_DIR)
        save_mols_to_sdf(
            unique_fragment_mols,
            output_file=os.path.join(
                settings.FRAGMENTS_OUTPUT_DIR, f"fragments.{target_id}.sdf"
            ),
        )
        log.info(f"Finished processing target: {target_id}")
    db_mgf_close_db_connection(db_mgf_connection)


if __name__ == "__main__":
    target_id_list = [
        # "CHEMBL1075189",
        # "CHEMBL1075126",
        # "CHEMBL3194",
        # "CHEMBL3864",
        # "CHEMBL3632452",
        # "CHEMBL288",
        # "CHEMBL5291543",
        "CHEMBL1741208",
    ]
    main(target_id_list)
