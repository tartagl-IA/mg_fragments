import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import scipy.cluster.hierarchy as sch
import sys
import streamlit as st
from rdkit.Chem import Draw, rdMolDescriptors
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity


ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(ROOT_DIR)

from db_mg_fragments.handlers import mols as db_mgf_mols_handler
from chem import filters as chem_filters
from chem import fragments as chem_fragments
from chem import utils as chem_utils

ss = st.session_state
if "reactive_patterns" not in ss:
    ss.reactive_patterns = chem_filters.get_reactive_patterns()
if "frag_mol_list_filtered" not in ss:
    ss.frag_mol_list_filtered = None

st.set_page_config(
    page_title="Molecule Explorer",
    page_icon="📊",
    layout="wide"
)

st.title("Molecule Explorer")

target_id_empty_value = "Select a Target ID"
target_id = st.sidebar.selectbox(
    label="Target ID",
    options=[target_id_empty_value] + db_mgf_mols_handler.get_available_targets()
)

with st.expander("Reactive Patterns"):
    data = pd.DataFrame([
        {
            "name": p["name"], 
            "smarts": f"`{p["smarts"]}`"
        }
        for p in ss.reactive_patterns
    ])
    data = data.astype({
        "name": "string",
        "smarts": "string"
    })
    
    st.table(data)
    
    index = st.number_input("index", min_value=0, step=1)
    name = st.text_input("name")
    smarts = st.text_input("smarts")
    button1, button2, button3, button4 = st.columns(4)
    with button1:
        if st.button("Add pattern"):
            try:
                ss.reactive_patterns.append({
                    "name": name,
                    "smarts": smarts,
                    "mol": chem_utils.mol_from_smarts(smarts)
                })
                st.rerun()
                st.info("Pattern added")
            except Exception as e:
                st.error(f"Error adding pattern: {e}")
    with button2:
        if st.button("Update pattern"):
            print(name, smarts)
            if index < len(ss.reactive_patterns):
                try:
                    ss.reactive_patterns[index] = {
                        "name": name if name else ss.reactive_patterns[index]["name"],
                        "smarts": smarts if smarts else ss.reactive_patterns[index]["smarts"],
                        "mol": chem_utils.mol_from_smarts(smarts if smarts else ss.reactive_patterns[index]["smarts"])
                    }
                    st.rerun()
                    st.info("Pattern updated")
                except Exception as e:
                    st.error(f"Error updating pattern: {e}")
    with button3:
        if st.button("Remove pattern"):
            try:
                if index < len(ss.reactive_patterns):
                    del ss.reactive_patterns[index]
                else:
                    st.error("Index out of range")
                st.rerun()
                st.info("Pattern removed")
            except Exception as e:
                st.error(f"Error removing pattern: {e}")
    with button4:
        if st.button("Save patterns"):
            try:
                reactive_patterns = [
                    {
                        "name": pattern["name"],
                        "smarts": pattern["smarts"]
                    }
                    for pattern in ss.reactive_patterns
                ]
                with open("chem/reactive_patterns.json", "w") as f:
                    import json
                    json.dump(reactive_patterns, f, indent=4)
                st.info("Patterns saved")
            except Exception as e:
                st.error(f"Error saving patterns: {e}")

if target_id is not target_id_empty_value:
    
    target_mols_data = [dict(data) for data in db_mgf_mols_handler.get_by_target(target_id)]
    for mol_data in target_mols_data:
        mol_data["mol"] = chem_utils.mol_from_smiles(mol_data["canonical_smiles"])
    
    total_molecule_num = st.sidebar.caption(f"Total molecule: {len(target_mols_data)}")
    
    st.sidebar.header("Molecule filters", divider=True)
    reactive = st.sidebar.toggle("Molecule reactive", value=False)

    filtered_mols_data = [
        mol_data
        for mol_data in target_mols_data
        if chem_filters.mol_reactive(mol_data["mol"], mol_patterns=[p["mol"] for p in ss.reactive_patterns]) is reactive
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
                                filters.append({
                                    "pattern": pattern["name"],
                                    "smarts": f"`{pattern["smarts"]}`",
                                    "found": res
                                })
                            st.table(filters)
                    st.divider()
                except Exception as e:
                    st.error(f"Error generating image for {mol_data['chembl_id']}: {e}")
                
            
            # img = Draw.MolsToGridImage(
            #     [mol_data["mol"] for mol_data in filtered_mols_data],
            #     legends=[f"{mol_data['canonical_smiles']} ({mol_data['chembl_id']})" for mol_data in filtered_mols_data],
            #     molsPerRow=4,
            #     subImgSize=(300, 300)
            # )
            # st.image(img)
    
    st.sidebar.header("Fragments settings", divider=True)
    fragment_min_dim = st.sidebar.number_input("Fragment minimum num atoms", min_value=1, step=1)
    fragment_max_dim = st.sidebar.number_input("Fragment maximum num atoms", min_value=0, step=1)
    fragment_flexibility = chem_filters.Flexibility(st.sidebar.pills("Fragment flexibility", chem_filters.Flexibility.values(), selection_mode="single", default=chem_filters.Flexibility.RIGID))
    fragment_max_num_rot_bonds = st.sidebar.number_input("Fragment max rotable bonds", min_value=1, step=1, disabled=fragment_flexibility!=chem_filters.Flexibility.FLEXIBLE)
    
    gen_fragment = st.sidebar.button("Generate fragments", icon="▶️")

    if gen_fragment:
        frag_list = set()
        for mol_data in filtered_mols_data:
            for bric in chem_fragments.brics_from_mol(mol_data["mol"]):
                 frag_list.add(bric)
        ss.frag_mol_list_filtered = [
            frag_mol
            for frag_mol in [chem_fragments.MolFromSmiles(frag) for frag in frag_list]
            if chem_filters.mol_dimension_range(frag_mol, fragment_min_dim, fragment_max_dim)
            if chem_filters.mol_flexibility(frag_mol, fragment_flexibility, fragment_max_num_rot_bonds)
        ]
        

if ss.frag_mol_list_filtered:
    st.write(f"Total fragments: {len(ss.frag_mol_list_filtered)}")
    
    threshold = st.slider(label="Cluster threshold", min_value=0.0, max_value=1.0, step=0.1)
    
    if st.button("Generate clusters", icon="▶️"):
    
        # Convert molecules to Morgan fingerprints
        fingerprints = [rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) for mol in ss.frag_mol_list_filtered]

        # Compute pairwise similarity matrix
        n_mols = len(ss.frag_mol_list_filtered)
        sim_matrix = np.zeros((n_mols, n_mols))

        for i in range(n_mols):
            for j in range(i+1, n_mols):
                sim_matrix[i, j] = TanimotoSimilarity(fingerprints[i], fingerprints[j])
                sim_matrix[j, i] = sim_matrix[i, j]  # Symmetric matrix

        # Convert similarity to distance (1 - similarity)
        dist_matrix = 1 - sim_matrix

        # Perform hierarchical clustering
        linkage_matrix = sch.linkage(dist_matrix, method="average")

        # Create clusters
        # threshold = 0.3  # Adjust clustering threshold
        clusters = sch.fcluster(linkage_matrix, threshold, criterion="distance")

        # Prepare clustering results
        clustered_mols = {}
        for idx, cluster_id in enumerate(clusters):
            if cluster_id not in clustered_mols:
                clustered_mols[cluster_id] = []
            clustered_mols[cluster_id].append(ss.frag_mol_list_filtered[idx])
        st.header(f"Clusters: {len(clustered_mols)}")
        
        # # Show Dendrogram
        # fig, ax = plt.subplots(figsize=(10, 5))
        # sch.dendrogram(linkage_matrix, labels=[f"Mol {i}" for i in range(n_mols)], leaf_rotation=90)
        # st.pyplot(fig)

        # # Show clustered molecules
        # for cluster_id, mols in clustered_mols.items():
        #     st.subheader(f"Cluster {cluster_id}")
        #     img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(30, 30))
        #     st.image(img)
        
        def mol_to_bytes(mol):
            from io import BytesIO
            img = Draw.MolToImage(mol, size=(300, 300))  # Generate molecule image
            buf = BytesIO()
            img.save(buf, format="PNG")
            return buf.getvalue()
        from rdkit.Chem import MolToSmiles
        for cluster_id, mols in clustered_mols.items():
            st.header(cluster_id)
            for mol in mols:
                st.image(mol_to_bytes(mol), caption=MolToSmiles(mol))
            break
        
        
