"""MG Fragments Streamlit App.

This module provides a Streamlit application for exploring and managing
molecular fragment databases, specifically the MGF and ChemDB databases.
It allows users to import molecular fragments from ChemDB to MGF and remove
fragments from the MGF database.
The application is designed to be user-friendly and provides a simple interface
for database operations.
"""

import logging

import streamlit as st

logging.getLogger("PIL.PngImagePlugin").setLevel(logging.WARNING)

# Set up Streamlit page configuration
st.set_page_config(page_title="Multi-page Streamlit App", page_icon="👋", layout="wide")

st.switch_page("pages/2_🔬_Molecule_Explorer.py")
