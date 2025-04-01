"""Importer module for importing mols from ChemDB to MGF DB."""

import sys

from db_chembl import close_db_connection as chembl_db_connection_close
from db_chembl import get_db_connection as chembl_db_connection_get
from db_chembl import utils as db_chembl_utils
from db_mg_fragments import close_db_connection as mgf_db_connection_close
from db_mg_fragments import get_db_connection as mgf_db_connection_get
from db_mg_fragments import handlers as mgf_db_apis
from logger import get_logger

log = get_logger(sys.modules[__name__].__name__)


def import_mol_by_targets_from_chem_db(target_id_list: list[str]) -> None:
    """Import mols associated to target_id_list from ChemDB to MGF DB.

    Args:
        target_id_list (list[str]): list of target IDs to import mols for
    """
    chembl_db_connection = chembl_db_connection_get()
    mgf_db_connection = mgf_db_connection_get()
    for i, target_id in enumerate(target_id_list, start=1):
        log.info(
            f"Importing mols associated to target {target_id} from ChemDB to MGF DB"
        )
        for mol in db_chembl_utils.get_mols_from_target_id(
            chembl_db_connection, target_id
        ):
            mgf_db_apis.mols.insert(mgf_db_connection, mol)
            log.info(
                f"Target ID: {target_id} [{i}/{len(target_id_list)}] \
                - Inserted mol: {mol['chembl_id']}"
            )
    chembl_db_connection_close(chembl_db_connection)
    mgf_db_connection_close(mgf_db_connection)
    log.info("All mols imported from ChemDB to MGF DB")


if __name__ == "__main__":
    import_mol_by_targets_from_chem_db(["CHEMBL1741208"])
