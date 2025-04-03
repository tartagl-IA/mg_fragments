import os
import sys

import pandas as pd
import streamlit as st
from rdkit.Chem import Draw

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(ROOT_DIR)

import settings
from chem import filters as chem_filters
from chem import fragments as chem_fragments
from chem import utils as chem_utils
from db_mg_fragments.handlers import mols as db_mgf_mols_handler

# --- setup ---

ss = st.session_state
if "reactive_patterns" not in ss:
    ss.reactive_patterns = chem_filters.get_reactive_patterns()
if "selected_target_id" not in ss:
    ss.selected_target_id = None
if "frag_mol_list_filtered" not in ss:
    ss.frag_mol_list_filtered = None
if "cluster_labels" not in ss:
    ss.cluster_labels = None
if "centroids" not in ss:
    ss.centroids = None

def reset():
    ss.frag_mol_list_filtered = None
    ss.cluster_labels = None
    ss.centroids = None


st.set_page_config(page_title="Molecule Explorer", page_icon="📊", layout="wide")


# --- page ---

st.title("Molecule Explorer")

target_id_empty_value = "Select a Target ID"
ss.selected_target_id = st.sidebar.selectbox(
    label="Target ID",
    options=[target_id_empty_value] + db_mgf_mols_handler.get_available_targets(),
    on_change=reset,
)

# Reactive patterns section

with st.expander("Reactive Patterns"):
    data = pd.DataFrame(
        [
            {"name": p["name"], "smarts": f"`{p["smarts"]}`"}
            for p in ss.reactive_patterns
        ]
    )
    data = data.astype({"name": "string", "smarts": "string"})

    st.table(data)

    st_index = st.number_input("index", min_value=0, step=1)
    st_name = st.text_input("name")
    st_smarts = st.text_input("smarts")
    st_button1, st_button2, st_button3, st_button4 = st.columns(4)
    button_state = None
    try:
        with st_button1:
            if st.button("Add pattern"):
                ss.reactive_patterns.append(
                    {
                        "name": st_name,
                        "smarts": st_smarts,
                        "mol": chem_utils.mol_from_smarts(st_smarts),
                    }
                )
                st.rerun()
                button_state = "Pattern added"
        with st_button2:
            if st.button("Update pattern"):
                if st_index < len(ss.reactive_patterns):
                    ss.reactive_patterns[st_index] = {
                        "name": (
                            st_name
                            if st_name
                            else ss.reactive_patterns[st_index]["name"]
                        ),
                        "smarts": (
                            st_smarts
                            if st_smarts
                            else ss.reactive_patterns[st_index]["smarts"]
                        ),
                        "mol": chem_utils.mol_from_smarts(
                            st_smarts
                            if st_smarts
                            else ss.reactive_patterns[st_index]["smarts"]
                        ),
                    }
                    st.rerun()
                    button_state = "Pattern updated"
        with st_button3:
            if st.button("Remove pattern"):
                if st_index < len(ss.reactive_patterns):
                    del ss.reactive_patterns[st_index]
                else:
                    st.error("Index out of range")
                st.rerun()
                button_state = "Pattern removed"
        with st_button4:
            if st.button("Save patterns"):
                reactive_patterns = [
                    {"name": pattern["name"], "smarts": pattern["smarts"]}
                    for pattern in ss.reactive_patterns
                ]
                with open("chem/reactive_patterns.json", "w") as f:
                    import json

                    json.dump(reactive_patterns, f, indent=4)
                button_state = "Patterns saved"
    except Exception as e:
        st.error(f"Error: {e}")
    if button_state:
        st.success(button_state)

st.divider()

# Fragments

if ss.selected_target_id is not target_id_empty_value:
    
    target_mols_data = [
        dict(data) for data in db_mgf_mols_handler.get_by_target(ss.selected_target_id)
    ]
    for mol_data in target_mols_data:
        mol_data["mol"] = chem_utils.mol_from_smiles(mol_data["canonical_smiles"])

    total_molecule_num = st.sidebar.caption(f"Total molecule: {len(target_mols_data)}")

    st.sidebar.header("Molecule filters", divider=True)
    reactive = st.sidebar.toggle("Molecule reactive", value=False)

    filtered_mols_data = [
        mol_data
        for mol_data in target_mols_data
        if chem_filters.mol_reactive(
            mol_data["mol"], mol_patterns=[p["mol"] for p in ss.reactive_patterns]
        )
        is reactive
    ]
    st.sidebar.caption(f"Filtered molecules: {len(filtered_mols_data)}")

    # remove counterions
    for mol_data in filtered_mols_data:
        max_mol = None
        max_mol_smiles = None
        for splitted_mol in mol_data["canonical_smiles"].split("."):
            if max_mol is None:
                max_mol = splitted_mol
                max_mol_smiles = splitted_mol
            else:
                if len(splitted_mol) > len(max_mol):
                    max_mol = splitted_mol
                    max_mol_smiles = splitted_mol
        if max_mol:
            mol_data["mol"] = chem_utils.mol_from_smiles(max_mol)
            mol_data["canonical_smiles"] = max_mol_smiles

    if filtered_mols_data:
        show = st.button("Show filtered molecules", icon="▶️")
        expose_patterns = st.toggle("Expose reactive patterns", value=False)
        
        if show:
            with st.expander("Filtered Molecules", expanded=True):
                st.subheader("Filtered molecules")
                for mol_data in filtered_mols_data:
                    try:
                        img_col, desc_col = st.columns([1, 1])
                        with img_col:
                            st.image(Draw.MolToImage(mol_data["mol"]))
                            st.write("CHEMBL ID: ", mol_data["chembl_id"])
                            st.write("SMILES: ", f'`{mol_data["canonical_smiles"]}`')
                        with desc_col:
                            if expose_patterns:
                                filters = []
                                for pattern in ss.reactive_patterns:
                                    if mol_data["mol"].HasSubstructMatch(pattern["mol"]):
                                        res = "✅"
                                    else:
                                        res = "❌"
                                    filters.append(
                                        {
                                            "pattern": pattern["name"],
                                            "smarts": f"`{pattern["smarts"]}`",
                                            "found": res,
                                        }
                                    )
                                st.table(filters)
                        st.divider()
                    except Exception as e:
                        st.error(f"Error generating image for {mol_data['chembl_id']}: {e}")

    st.sidebar.header("Fragments settings", divider=True)
    fragment_min_dim = st.sidebar.number_input(
        "Fragment minimum num atoms", min_value=1, step=1, value=12
    )
    fragment_max_dim = st.sidebar.number_input(
        "Fragment maximum num atoms", min_value=0, step=1
    )
    fragment_flexibility = chem_filters.Flexibility(
        st.sidebar.pills(
            "Fragment flexibility",
            chem_filters.Flexibility.values(),
            selection_mode="single",
            default=chem_filters.Flexibility.RIGID,
        )
    )
    fragment_max_num_rot_bonds = st.sidebar.number_input(
        "Fragment max rotable bonds",
        min_value=1,
        step=1,
        disabled=fragment_flexibility != chem_filters.Flexibility.FLEXIBLE,
    )

    if st.sidebar.button("Generate fragments", icon="▶️"):
        with st.spinner("Generating fragments..."):
            frag_list = set()
            for mol_data in filtered_mols_data:
                for bric in chem_fragments.brics_from_mol(mol_data["mol"]):
                    frag_list.add(bric)
            ss.frag_mol_list_filtered = [
                frag_mol
                for frag_mol in [chem_fragments.MolFromSmiles(frag) for frag in frag_list]
                if chem_filters.mol_dimension_range(
                    frag_mol, fragment_min_dim, fragment_max_dim
                )
                if chem_filters.mol_flexibility(
                    frag_mol, fragment_flexibility, fragment_max_num_rot_bonds
                )
            ]

    if ss.frag_mol_list_filtered:
        st.write(f"Total fragments: {len(ss.frag_mol_list_filtered)}")
        if st.button("Show fragments", icon="🙈"):
            with st.expander("Fragments"):
                for frag_mol in ss.frag_mol_list_filtered:
                    try:
                        st.image(Draw.MolToImage(frag_mol))
                    except Exception as e:
                        st.error(f"Error generating image: {e}")


if ss.frag_mol_list_filtered:
    
    st.divider()
    
    clustering_type = st.selectbox(
        "Clustering type",
        ["", "Maximum Common Substructure", "Tanimoto Similarity"],
    )
    cluster_method = st.selectbox("Cluster method", ["", "distance", "maxclust"])
    if cluster_method == "distance":
        t = st.number_input(
            "Cluster threshold", min_value=0.0, max_value=5.0, step=0.1, value=0.3
        )
    elif cluster_method == "maxclust":
        t = st.number_input("Maximum number of clusters", min_value=1, step=1)
    else:
        t = None

    if clustering_type and cluster_method:

        clustering_method = None
        if clustering_type == "Maximum Common Substructure":
            from chem.clustering import mcs as clustering_method
        elif clustering_type == "Tanimoto Similarity":
            from chem.clustering import tanimoto as clustering_method

        if clustering_method:
            if st.button("Generate clusters", icon="▶️"):
                with st.spinner("Generating clusters..."):
                    ss.cluster_labels = clustering_method.hierarchical_clustering(
                        ss.frag_mol_list_filtered, cluster_method, t
                    )
                    ss.centroids = clustering_method.find_cluster_centroids(
                        ss.frag_mol_list_filtered, ss.cluster_labels
                    )
                    st.toast("Prepared clustering results", icon="🎉")

                unique_cluster_labels = set(ss.cluster_labels)
                st.subheader(
                    f"Number of generated clusters: {len(unique_cluster_labels)}"
                )

        if ss.cluster_labels is not None and ss.centroids is not None:

            if st.button("Save to SDF file", icon="⬇️"):
                with st.spinner("Saving to SDF file..."):
                    from chem.utils import save_mols_to_sdf

                    save_mols_to_sdf(
                        [
                            ss.centroids[cluster_id]
                            for cluster_id in set(ss.cluster_labels)
                        ],
                        output_file=os.path.join(
                            ROOT_DIR,
                            settings.FRAGMENTS_OUTPUT_DIR,
                            f"{ss.selected_target_id}_centroids_fragments.sdf",
                        ),
                    )
                st.toast("Saved to SDF file", icon="🎉")

            if st.button("Show clustered molecules", icon="🔍"):

                for current_cluster_id in set(ss.cluster_labels):
                    st.subheader(f"Cluster {current_cluster_id}")
                    st.image(
                        chem_utils.mol_to_bytes(ss.centroids[current_cluster_id]),
                        caption=chem_utils.smiles_from_mol(
                            ss.centroids[current_cluster_id]
                        ),
                    )
                    idx_mol_in_cluster = [
                        idx
                        for idx, cluster_id in enumerate(ss.cluster_labels)
                        if cluster_id == current_cluster_id
                    ]
                    with st.expander(
                        f"{len(idx_mol_in_cluster)} molecule in cluster", expanded=False
                    ):
                        for idx in idx_mol_in_cluster:
                            mol = ss.frag_mol_list_filtered[idx]
                            st.image(
                                chem_utils.mol_to_bytes(mol),
                                caption=chem_utils.smiles_from_mol(mol),
                            )
                    st.divider()
