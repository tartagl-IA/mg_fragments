"""Database Explorer for MGF and ChemDB.

This module provides a Streamlit application for exploring and managing
molecular fragment databases, specifically the MGF and ChemDB databases.
It allows users to import molecular fragments from ChemDB to MGF and remove
fragments from the MGF database.
The application is designed to be user-friendly and provides a simple interface
for database operations.
"""

import os
import sys

import streamlit as st

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(ROOT_DIR)

import db_chembl
import db_chembl.utils as db_chembl_utils
import db_mg_fragments
import db_mg_fragments.handlers.mols as db_mgf_mols_handlers

ss = st.session_state
sb = st.sidebar
if "db_chembl_target_id_list" not in ss:
    ss.db_chembl_target_id_list = None
if "db_mgf_target_id_list" not in ss:
    ss.db_mgf_target_id_list = db_mgf_mols_handlers.get_available_targets()


def import_mol_by_targets_from_chem_db(target_id_list: list[str]) -> None:
    """Import mols associated to target_id_list from ChemDB to MGF DB.

    Args:
        target_id_list (list[str]): list of target IDs to import mols for
    """
    db_chembl_connection = db_chembl.get_db_connection()
    mgf_db_connection = db_mg_fragments.get_db_connection()
    for i, target_id in enumerate(target_id_list, start=1):
        st.toast(
            f"Importing mols associated to target {target_id} from ChemDB to MGF DB"
        )
        for mol in db_chembl_utils.get_mols_from_target_id(
            db_chembl_connection, target_id
        ):
            db_mgf_mols_handlers.insert(mgf_db_connection, mol)
            st.toast(
                f"Target ID: {target_id} [{i}/{len(target_id_list)}] \
                - Inserted mol: {mol['chembl_id']}"
            )
    db_chembl_connection.close()
    mgf_db_connection.close()
    st.info("All mols imported from ChemDB to MGF DB")


st.set_page_config(page_title="DB Explorer", page_icon="📊", layout="wide")
st.title("DB Explorer")

action = sb.selectbox("Actions", ["Import target", "Remove Target"])
with st.expander("DB MGF Available Target IDs", expanded=True):
    st.table(ss.db_mgf_target_id_list)

if action == "Import target":
    st.subheader("Import target from CHEMBL DB")

    if st.button("Import target from CHEMBL DB", icon="🔄"):
        with st.spinner("Retrieving Target"):
            try:
                from db_chembl.utils import get_available_target_ids

                ss.db_chembl_target_id_list = get_available_target_ids()
                st.success("Target IDs retrieved successfully.")
            except Exception as e:
                st.error(f"Error retrieving target IDs: {e}")

    if ss.db_chembl_target_id_list:
        selected_target_id = st.selectbox(
            "Available Target IDs", [""] + ss.db_chembl_target_id_list
        )
        if selected_target_id != "":
            if st.button(f"Import Target {selected_target_id}"):
                with st.spinner("Importing Target"):
                    try:
                        import_mol_by_targets_from_chem_db([selected_target_id])
                        st.success(
                            f"Target ID {selected_target_id} imported successfully."
                        )
                    except Exception as e:
                        st.error(f"Error importing target ID {selected_target_id}: {e}")

if ss.db_mgf_target_id_list and action == "Remove Target":
    st.subheader("Remove Target from DB MGF")

    selected_target_id = st.selectbox(
        "Available Target IDs", [""] + ss.db_mgf_target_id_list
    )
    if selected_target_id:
        if st.button(f"Remove Target {selected_target_id}"):
            with st.spinner("Removing Target"):
                try:
                    connection = db_mg_fragments.get_db_connection()
                    db_mgf_mols_handlers.remove_by_target_id(
                        connection, selected_target_id
                    )
                    connection.close()
                    st.success(f"Target ID {selected_target_id} removed successfully.")
                except Exception as e:
                    st.error(f"Error removing target ID {selected_target_id}: {e}")
