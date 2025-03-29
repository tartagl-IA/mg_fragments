"""Handler for the 'mols' table in the SQLite database."""

import sqlite3
from typing import Any, Generator

from logger import get_logger

TABLE_NAME = "mols"
log = get_logger("DB MGF")


def insert(connection: sqlite3.Connection, mol: dict[str, Any]) -> None:
    """Inserts a molecule into the 'mols' table.

    Args:
        connection (sqlite3.Connection): SQLite database connection.
        mol (dict): Dictionary containing molecule data.

    Returns:
        None
    """
    log.debug(f"Inserting into '{TABLE_NAME}' table")
    cursor = connection.cursor()
    query = f"""
        INSERT INTO {TABLE_NAME} (
            target_id,
            chembl_id,
            canonical_smiles
        ) VALUES (?, ?, ?)
    """
    data = (
        mol["target_id"],
        mol["chembl_id"],
        mol["canonical_smiles"],
    )
    cursor.execute(query, data)
    connection.commit()
    log.debug(f"Inserted into '{TABLE_NAME}' table")


def get_by_target(
    connection: sqlite3.Connection, target_id: str
) -> Generator[dict[str, Any], None, None]:
    """Retrieves a molecule from the 'mols' table based on target_id.

    Args:
        connection (sqlite3.Connection): SQLite database connection.
        target_id (str): Target ID of the molecule.

    Returns:
        dict: Dictionary containing molecule data.
    """
    log.debug(f"Fetching from '{TABLE_NAME}' table by target_id: {target_id}")
    cursor = connection.cursor()
    query = f"""
        SELECT * FROM {TABLE_NAME} WHERE target_id = ?
    """
    cursor.execute(query, (target_id,))
    for row in cursor:
        yield row
    log.debug(f"Fetched all result for target ID: {target_id}")
