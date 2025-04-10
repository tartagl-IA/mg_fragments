"""Molecule Explorer.

This is a Streamlit app for exploring molecules and their fragments.
It allows users to filter molecules based on reactivity patterns,
generate fragments, and perform clustering on the fragments.
"""

import json
import os
import sys
from dataclasses import dataclass

import pandas as pd
import streamlit as st
from rdkit.Chem import Draw, Mol

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(ROOT_DIR)

import db_mg_fragments.handlers.mols as db_mgf_mols_handler
import settings
from chem import clustering as chem_clustering
from chem import filters as chem_filters
from chem import fragments as chem_fragments
from chem import utils as chem_utils


@dataclass
class MoleculeData(db_mgf_mols_handler.Molecule):
    """Molecule data class."""

    mol: Mol


# --- setup ---

ss = st.session_state
if "selected_target_id" not in ss:
    ss.selected_target_id = None
if "target_mols_data" not in ss:
    ss.target_mols_data = None
if "target_mols_data_filtered" not in ss:
    ss.target_mols_data_filtered = None
if "reactive_toggle" not in ss:
    ss.reactive_toggle = None
if "reactive_pattern_list" not in ss:
    ss.reactive_pattern_list = None
if "db_mgf_target_id_list" not in ss:
    ss.db_mgf_target_id_list = None
if "selected_target_id_idx" not in ss:
    ss.selected_target_id_idx = 0
if "fragment_flexibility" not in ss:
    ss.fragment_flexibility = None
if "frag_mol_list_filtered" not in ss:
    ss.frag_mol_list_filtered = None
if "cluster_labels" not in ss:
    ss.cluster_labels = None
if "centroids" not in ss:
    ss.centroids = None


def selected_target_id_on_change():
    """On change of the selected target ID, reset the filtered molecules and fragments."""
    ss.target_mols_data = None
    ss.target_mols_data_filtered = None
    ss.frag_mol_list_filtered = None
    ss.cluster_labels = None
    ss.centroids = None


def reactive_toggle_on_change():
    """On change of the reactive toggle, reset the filtered molecules and fragments."""
    ss.target_mols_data_filtered = None
    ss.frag_mol_list_filtered = None
    ss.cluster_labels = None
    ss.centroids = None

    if ss.reactive_toggle:
        st.toast("Filtered reactive molecules", icon="🔄")
    elif ss.reactive_toggle is False:
        st.toast("Filtered not-reactive molecules", icon="🔄")


# --- page ---


st.set_page_config(page_title="Molecule Explorer", page_icon="📊", layout="wide")

with st.popover("Session State", icon="🔄"):
    st.write(ss.cluster_labels)
    st.write(ss.centroids)

if ss.reactive_pattern_list is None:
    ss.reactive_pattern_list = chem_filters.get_reactive_pattern_list()
if ss.db_mgf_target_id_list is None:
    ss.db_mgf_target_id_list = db_mgf_mols_handler.get_available_targets()

st.title("Molecule Explorer")

ss.selected_target_id = st.sidebar.selectbox(
    label="Target ID",
    options=ss.db_mgf_target_id_list,
    on_change=selected_target_id_on_change,
    # index=ss.selected_target_id_idx
)

if ss.selected_target_id:
    # st.toast(f"Selected target ID: {ss.selected_target_id}", icon="🔄")
    # if ss.target_mols_data is None:
    ss.target_mols_data = []
    for data in db_mgf_mols_handler.get_by_target(ss.selected_target_id):
        ss.target_mols_data.append(
            MoleculeData(
                target_id=data["target_id"],
                chembl_id=data["chembl_id"],
                canonical_smiles=data["canonical_smiles"],
                mol=chem_utils.mol_from_smiles(data["canonical_smiles"]),
            )
        )

    total_molecule_num = st.sidebar.write(
        f"Total molecule: `{len(ss.target_mols_data)}`"
    )

    # Reactive patterns section

    st.subheader("1. Reactive filter", divider=True)
    with st.expander("Reactive Patterns", icon="🧩"):
        data = pd.DataFrame(
            [
                {"name": p.name, "smarts": f"`{p.smarts}`"}
                for p in ss.reactive_pattern_list
            ]
        )

        st.table(data)

        st_index = st.number_input("index", min_value=0, step=1)
        st_name = st.text_input("name")
        st_smarts = st.text_input("smarts")
        st_button1, st_button2, st_button3, st_button4 = st.columns(4)
        button_state = None
        try:
            with st_button1:
                if st.button("Add pattern"):
                    ss.reactive_pattern_list.append(
                        chem_filters.ReactivePattern(
                            name=st_name,
                            smarts=st_smarts,
                            mol=chem_utils.mol_from_smarts(st_smarts),
                        )
                    )
                    st.rerun()
                    button_state = "Pattern added"
            with st_button2:
                if st.button("Update pattern"):
                    if st_index < len(ss.reactive_pattern_list):
                        if st_name:
                            ss.reactive_pattern_list[st_index].name = st_name
                        if st_smarts:
                            ss.reactive_pattern_list[st_index].smarts = st_smarts
                            ss.reactive_pattern_list[st_index].mol = (
                                chem_utils.mol_from_smarts(st_smarts)
                            )
                        st.rerun()
                        button_state = "Pattern updated"
            with st_button3:
                if st.button("Remove pattern"):
                    if st_index < len(ss.reactive_pattern_list):
                        del ss.reactive_pattern_list[st_index]
                    else:
                        st.error("Index out of range")
                    st.rerun()
                    button_state = "Pattern removed"
            with st_button4:
                if st.button("Save patterns"):
                    with open("chem/reactive_patterns.json", "w") as f:
                        json.dump(
                            [rp.to_dict() for rp in ss.reactive_pattern_list],
                            f,
                            indent=4,
                        )
                    button_state = "Patterns saved"
        except Exception as e:
            st.error(f"Error: {e}")
        if button_state:
            st.success(button_state)

if ss.target_mols_data:
    ss.reactive_toggle = st.toggle(
        "Molecule reactive", value=False, on_change=reactive_toggle_on_change
    )
    ss.target_mols_data_filtered = [
        mol_data
        for mol_data in ss.target_mols_data
        if chem_filters.mol_reactive(
            mol_data.mol, reactive_pattern_list=ss.reactive_pattern_list
        )
        is ss.reactive_toggle
    ]
    # remove counterions
    for mol_data in ss.target_mols_data_filtered:
        mol_data.canonical_smiles = chem_utils.remove_counterions_from_smiles(
            mol_data.canonical_smiles
        )
        mol_data.mol = chem_utils.mol_from_smiles(mol_data.canonical_smiles)
    st.sidebar.write(f"Filtered molecules: `{len(ss.target_mols_data_filtered)}`")
    c1, c2, _ = st.columns([1, 1, 4])
    with c1:
        show_mol_filtered = st.button("Show filtered molecules", icon="▶️")
    with c2:
        expose_patterns = st.toggle(
            "Expose reactive patterns", value=ss.reactive_toggle
        )
    if show_mol_filtered:
        with st.spinner("Generating filtered molecules..."):
            with st.expander("Filtered Molecules", expanded=True):
                st.subheader("Filtered molecules")
                for mol_data in ss.target_mols_data_filtered:
                    try:
                        img_col, desc_col = st.columns([1, 1])
                        with img_col:
                            st.image(Draw.MolToImage(mol_data.mol))
                            st.write("CHEMBL ID: ", mol_data.chembl_id)
                            st.write("SMILES: ", f"`{mol_data.canonical_smiles}`")
                        with desc_col:
                            if expose_patterns:
                                filters = []
                                for pattern in ss.reactive_pattern_list:
                                    if mol_data.mol.HasSubstructMatch(pattern.mol):
                                        res = "✅"
                                    else:
                                        res = "❌"
                                    filters.append(
                                        {
                                            "pattern": pattern.name,
                                            "smarts": f"`{pattern.smarts}`",
                                            "found": res,
                                        }
                                    )
                                st.table(filters)
                        st.divider()
                    except Exception as e:
                        st.error(
                            f"Error generating image for {mol_data.chembl_id}: {e}"
                        )

# Fragments

if ss.target_mols_data_filtered:
    st.subheader("2. Fragments generation", divider=True)

    with st.expander("Fragment generation settings", expanded=False):
        fragment_min_dim = st.number_input(
            "Fragment minimum num atoms", min_value=1, step=1, value=12
        )
        fragment_max_dim = st.number_input(
            "Fragment maximum num atoms", min_value=0, step=1, value=60
        )
        ss.fragment_flexibility = chem_filters.Flexibility(
            st.pills(
                "Fragment flexibility",
                chem_filters.Flexibility.values(),
                selection_mode="single",
                default=chem_filters.Flexibility.RIGID,
            )
        )
        fragment_max_num_rot_bonds = st.number_input(
            "Fragment max rotable bonds",
            min_value=1,
            step=1,
            disabled=ss.fragment_flexibility != chem_filters.Flexibility.FLEXIBLE,
        )

    c1, c2, _ = st.columns([1, 1, 4])
    show_fragments = False
    with c1:
        if st.button("Generate fragments", icon="▶️"):
            with st.spinner("Generating fragments..."):
                frag_list = set()
                for mol_data in ss.target_mols_data_filtered:
                    for bric in chem_fragments.brics_from_mol(mol_data.mol):
                        frag_list.add(bric)
                ss.frag_mol_list_filtered = [
                    frag_mol
                    for frag_mol in [
                        chem_fragments.MolFromSmiles(frag) for frag in frag_list
                    ]
                    if chem_filters.mol_dimension_range(
                        frag_mol, fragment_min_dim, fragment_max_dim
                    )
                    if chem_filters.mol_flexibility(
                        frag_mol, ss.fragment_flexibility, fragment_max_num_rot_bonds
                    )
                ]
                st.toast("Fragments generated", icon="🎉")
    with c2:
        if ss.frag_mol_list_filtered:
            show_fragments = st.button("Show fragments", icon="🔍")

    if ss.frag_mol_list_filtered:
        st.sidebar.write(f"Generated fragments: `{len(ss.frag_mol_list_filtered)}`")

    if ss.frag_mol_list_filtered and show_fragments:
        with st.spinner("Generating fragments images..."):
            with st.expander("Generated fragments images", expanded=True):
                for frag_mol in ss.frag_mol_list_filtered:
                    try:
                        st.image(Draw.MolToImage(frag_mol))
                    except Exception as e:
                        st.error(f"Error generating image: {e}")


if ss.frag_mol_list_filtered:

    st.subheader("3. Clustering", divider=True)

    with st.expander("Clustering settings", expanded=True):
        clustering_type_module = None
        clustering_type = st.selectbox(
            "Clustering type",
            [clust_type.value for clust_type in chem_clustering.ClusteringType],
            index=None,
            placeholder="Select clustering type",
        )
        if clustering_type == chem_clustering.ClusteringType.MCS.value:
            from chem.clustering import mcs as clustering_type_module
        elif clustering_type == chem_clustering.ClusteringType.TANIMOTO.value:
            from chem.clustering import tanimoto as clustering_type_module

        cluster_method = st.selectbox(
            "Cluster method",
            [clust_method.value for clust_method in chem_clustering.ClusteringMethod],
            index=None,
            placeholder="Select cluster method",
        )
        if cluster_method == chem_clustering.ClusteringMethod.DIST.value:
            t = st.number_input(
                "Cluster threshold", min_value=0.0, max_value=5.0, step=0.1, value=0.3
            )
        elif cluster_method == chem_clustering.ClusteringMethod.MAX_CLUSTERS.value:
            t = st.number_input(
                "Maximum number of clusters", min_value=1, step=1, value=6
            )
        else:
            t = None

        st.toggle("Sanitize molecules", value=False, key="sanitize_molecules")

    if clustering_type_module and cluster_method:

        if st.button("Generate clusters", icon="▶️"):
            if ss.sanitize_molecules:
                frag_mol_list_filtered = [
                    chem_utils.sanitise_mol(mol) for mol in ss.frag_mol_list_filtered
                ]
                st.toast("Sanitized molecules", icon="🎉")
            else:
                frag_mol_list_filtered = ss.frag_mol_list_filtered
            with st.spinner("Generating clusters..."):
                for item in clustering_type_module.hierarchical_clustering(
                    frag_mol_list_filtered, cluster_method, t
                ):
                    if "log" in item:
                        st.toast(item["log"])
                    elif "result" in item:
                        ss.cluster_labels = item["result"]
                for item in clustering_type_module.find_cluster_centroids(
                    frag_mol_list_filtered, ss.cluster_labels
                ):
                    if "log" in item:
                        st.toast(item["log"])
                    elif "result" in item:
                        ss.centroids = item["result"]
                st.toast("Prepared clustering results", icon="🎉")
    if ss.cluster_labels is not None:
        st.toast(
            f"Generated {len(set(ss.cluster_labels))} total clusters",
            icon="🎉",
        )
        st.sidebar.write(
            f"Number of generated clusters: `{len(set(ss.cluster_labels))}`"
        )

    if ss.cluster_labels is not None:  # and ss.centroids is not None:

        c1, c2, _ = st.columns([1, 1, 4])
        show_clusters = None
        with c1:
            show_clusters = st.button("Show clustered molecules", icon="🔍")
        with c2:
            with st.popover("Save to SDF", icon="💾", use_container_width=True):
                st.subheader("Save clustered molecules to SDF file")
                output_dir = os.path.join(
                    ROOT_DIR,
                    settings.FRAGMENTS_OUTPUT_DIR,
                )
                st.write(f"Output directory: `{output_dir}`")
                st.markdown(
                    f"""
                Available keywords:
                |Keyword|Current value|
                |---|---|
                |`target_id`|{ss.selected_target_id}|
                |`reactive`|{ss.reactive_toggle}|
                |`flexibility`|{chem_filters.Flexibility(ss.fragment_flexibility).value}|
                """
                )
                st.info(
                    "Use the keywords to create a custom file name with the format: `{target_id}_centroids_fragments.sdf`"
                )
                output_file_name_preformatted = "{target_id}_{reactive}_{flexibility}_centroids_fragments.sdf".format(
                    target_id=ss.selected_target_id,
                    reactive="reactive" if ss.reactive_toggle else "non-reactive",
                    flexibility=chem_filters.Flexibility(ss.fragment_flexibility).value,
                )
                output_file_name = st.text_input(
                    "Output file name",
                    value=output_file_name_preformatted,
                )
                output_file_name = output_file_name.format(
                    target_id=ss.selected_target_id,
                    reactive="reactive" if ss.reactive_toggle else "non-reactive",
                    flexibility=chem_filters.Flexibility(ss.fragment_flexibility).value,
                )
                if not output_file_name.endswith(".sdf"):
                    output_file_name += ".sdf"
                if st.button("Save to SDF file", icon="⬇️"):
                    chem_utils.save_mols_to_sdf(
                        [
                            ss.centroids[cluster_id]
                            for cluster_id in set(ss.cluster_labels)
                        ],
                        output_file=os.path.join(output_dir, output_file_name),
                    )
                    st.info("Saved to SDF file", icon="🎉")

        if show_clusters:
            with st.spinner("Generating clustered molecules..."):
                for current_cluster_id in set(ss.cluster_labels):
                    st.subheader(f"Cluster {current_cluster_id}")
                    if ss.centroids and current_cluster_id in ss.centroids:
                        st.image(
                            chem_utils.mol_to_bytes(ss.centroids[current_cluster_id]),
                            caption=chem_utils.smiles_from_mol(
                                ss.centroids[current_cluster_id]
                            ),
                        )
                    else:
                        st.write("No centroid found for this cluster")
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
