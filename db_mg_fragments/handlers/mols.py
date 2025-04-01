"""Handler for the 'mols' table in the SQLite database."""

import sqlite3
from typing import Any, Generator

from logger import get_logger

from .. import get_db_connection

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


def get_available_targets() -> list[str]:
    log.debug(f"Fetching from '{TABLE_NAME}' table available targets")
    connection = get_db_connection()
    cursor = connection.cursor()
    query = f"""
        SELECT DISTINCT target_id FROM mols
    """
    cursor.execute(query)
    res = [row["target_id"] for row in cursor.fetchall()]
    connection.close()
    return res


def get_by_target(target_id: str) -> Generator[dict[str, Any], None, None]:
    """Retrieves a molecule from the 'mols' table based on target_id.

    Args:
        connection (sqlite3.Connection): SQLite database connection.
        target_id (str): Target ID of the molecule.

    Returns:
        dict: Dictionary containing molecule data.
    """
    log.debug(f"Fetching from '{TABLE_NAME}' table by target_id: {target_id}")
    connection = get_db_connection()
    cursor = connection.cursor()
    query = f"""
        SELECT * FROM {TABLE_NAME} WHERE target_id = ?
    """
    cursor.execute(query, (target_id,))
    for row in cursor:
        yield row
    connection.close()
    log.debug(f"Fetched all result for target ID: {target_id}")
