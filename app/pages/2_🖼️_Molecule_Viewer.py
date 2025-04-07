"""Molecule Viewer.

This script allows users to view molecules from a selected SDF file.
It uses RDKit to read the SDF file and Streamlit to display the molecules.
"""

import os
import sys
from io import BytesIO

import streamlit as st
from rdkit.Chem import Draw, Mol, MolToSmiles, SDMolSupplier

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(ROOT_DIR)

from settings import FRAGMENTS_OUTPUT_DIR


def get_available_targets() -> list[str]:
    """Retrieve the list of available target files.

    Returns:
        list[str]: available SDF target files
    """
    return [
        sdf_file
        for sdf_file in os.listdir(os.path.join(os.getcwd(), FRAGMENTS_OUTPUT_DIR))
        if sdf_file.endswith(".sdf")
    ]


def get_images(sdf_file) -> list[Mol]:
    """Get images of molecules from the specified SDF file.

    Args:
        sdf_file (str): The name of the SDF file.

    Returns:
        list[Mol]: List of valid RDKit molecules.
    """
    supplier = SDMolSupplier(os.path.join(ROOT_DIR, FRAGMENTS_OUTPUT_DIR, sdf_file))
    return [mol for mol in supplier if mol is not None]


def mol_to_bytes(mol: Mol) -> bytes:
    """Convert RDKit molecule to bytes for display in Streamlit.

    Args:
        mol (Mol): RDKit molecule object.

    Returns:
        bytes: Image bytes of the molecule.
    """
    img = Draw.MolToImage(mol, size=(300, 300))
    buf = BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


st.set_page_config(page_title="Molecule Viewer", page_icon="🧪", layout="wide")
st.title("Molecule Viewer")

file = st.selectbox(label="Target File", options=get_available_targets())
if file:
    molecules = get_images(file)
    if molecules:
        for mol in molecules:
            st.image(mol_to_bytes(mol), caption=MolToSmiles(mol))
    else:
        st.warning("No valid molecules found in the file.")
