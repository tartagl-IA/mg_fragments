import os
import sys
from io import BytesIO

import streamlit as st
from rdkit.Chem import Draw, MolToSmiles, SDMolSupplier

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(ROOT_DIR)

from settings import FRAGMENTS_OUTPUT_DIR


def get_available_targets() -> list[str]:
    return os.listdir(os.path.join(os.getcwd(), FRAGMENTS_OUTPUT_DIR))


def get_images(sdf_file):
    supplier = SDMolSupplier(os.path.join(ROOT_DIR, FRAGMENTS_OUTPUT_DIR, sdf_file))

    # Filter valid molecules
    return [mol for mol in supplier if mol is not None]


# Function to convert RDKit molecule to an image
def mol_to_bytes(mol):
    img = Draw.MolToImage(mol, size=(300, 300))  # Generate molecule image
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
