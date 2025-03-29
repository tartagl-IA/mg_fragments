"""Database connection module for MG Fragment database."""

import os
import sqlite3

DB_NAME = "db_mg_fragments.db"
DB_PATH = os.path.join(os.path.dirname(__file__), DB_NAME)


def get_db_connection() -> sqlite3.Connection:
    """Get a connection to the SQLite database.

    Returns:
        sqlite3.Connection: A connection object to the SQLite database.
    """
    connection = sqlite3.connect(DB_PATH)
    connection.row_factory = sqlite3.Row
    return connection


def close_db_connection(connection: sqlite3.Connection) -> None:
    """Close the SQLite database connection.

    Args:
        connection (sqlite3.Connection): A connection object to the SQLite database.
    """
    connection.close()
