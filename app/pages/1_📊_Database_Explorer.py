import os
import sys

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(ROOT_DIR)

import streamlit as st

ss = st.session_state
sb = st.sidebar
if "db_chembl_target_id_list" not in ss:
    ss.db_chembl_target_id_list = None
if "db_mgf_target_id_list" not in ss:
    from db_mg_fragments.handlers.mols import get_available_targets
    ss.db_mgf_target_id_list = get_available_targets()

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
        if selected_target_id:
            if st.button(f"Import Target {selected_target_id}"):
                with st.spinner("Importing Target"):
                    try:
                        from importer import import_mol_by_targets_from_chem_db
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
                    from db_mg_fragments import get_db_connection
                    from db_mg_fragments.handlers.mols import remove_by_target_id
                    connection = get_db_connection()
                    remove_by_target_id(connection, selected_target_id)
                    connection.close()
                    st.success(
                        f"Target ID {selected_target_id} removed successfully."
                    )
                except Exception as e:
                    st.error(f"Error removing target ID {selected_target_id}: {e}")
    