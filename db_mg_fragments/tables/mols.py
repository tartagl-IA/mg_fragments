"""Module to create and manage the 'mols' table in the database."""

import sqlite3

from logger import get_logger

log = get_logger("DB MGF")

TABLE_NAME = "mols"


def create(connection: sqlite3.Connection) -> None:
    """Creates the 'mols' table in the database.

    Args:
        connection (sqlite3.Connection): The SQLite connection object.
    """
    log.debug(f"Creating '{TABLE_NAME}' table")
    cursor = connection.cursor()
    query = f"""
        CREATE TABLE IF NOT EXISTS {TABLE_NAME} (
            target_id TEXT,
            chembl_id TEXT,
            canonical_smiles TEXT,
            PRIMARY KEY (target_id, chembl_id)
        )
    """
    cursor.execute(query)
    connection.commit()
    log.debug(f"Table '{TABLE_NAME}' created")
