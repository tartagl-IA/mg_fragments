"""Logging module for the project."""

import logging

# Configure logging
logging.basicConfig(
    format="[%(asctime)s][%(levelname)-7s][%(module)s] %(message)s",
    level=logging.DEBUG,
    datefmt="%Y-%m-%d %H:%M:%S",
)


def get_logger(name: str) -> logging.Logger:
    """Get a logger with the specified name.

    Args:
        name (str): logger name

    Returns:
        logging.Logger: logger instance
    """
    return logging.getLogger(name)
