"""Utilities for the db_chembl package."""

import sqlite3
from typing import Any, Generator

from logger import get_logger

log = get_logger("DB CHEMBL")


def get_mols_from_target_id(
    connection: sqlite3.Connection, target_id: str
) -> Generator[dict[str, Any], None, None]:
    """Fetches all molecules from the 'mols' table for a given target ID.

    Args:
        connection (sqlite3.Connection): SQLite database connection.
        target_id (str): The target ID to filter molecules.

    Returns:
        list: List of tuples containing molecule data.
    """
    log.debug(f"Fetching molecules for target ID: {target_id}")
    cursor = connection.cursor()
    query = """
        SELECT DISTINCT
            td.chembl_id AS target_id,
            md.chembl_id,
            cs.canonical_smiles
        FROM target_dictionary td
        JOIN assays ass ON ass.tid = td.tid
        JOIN activities act ON act.assay_id = ass.assay_id
        JOIN molecule_dictionary md ON md.molregno = act.molregno
        JOIN compound_structures cs ON cs.molregno = act.molregno
        WHERE td.chembl_id = ?
    """
    cursor.execute(query, (target_id,))
    for row in cursor:
        yield row
    log.debug(f"Fetched all result for target ID: {target_id}")
