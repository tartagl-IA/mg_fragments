"""MG Fragments Streamlit App.

This module provides a Streamlit application for exploring and managing
molecular fragment databases, specifically the MGF and ChemDB databases.
It allows users to import molecular fragments from ChemDB to MGF and remove
fragments from the MGF database.
The application is designed to be user-friendly and provides a simple interface
for database operations.
"""

import logging
import os
import sys

import streamlit as st

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(ROOT_DIR)

logging.getLogger("PIL.PngImagePlugin").setLevel(logging.WARNING)

# Set up Streamlit page configuration
st.set_page_config(page_title="MG Fragments", page_icon="🧪", layout="wide")

with open("README.md", "r") as file:
    readme = file.read()
st.markdown(readme, unsafe_allow_html=True)
