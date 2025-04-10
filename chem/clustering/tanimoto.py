"""This module provides functions for hierarchical clustering of RDKit molecules using Tanimoto similarity.

It includes functions to generate Morgan fingerprints, compute pairwise Tanimoto similarity,
and find cluster centroids based on the computed clusters.
The clustering is performed using SciPy's hierarchical clustering methods.
The module is designed to work with molecular structures represented as RDKit Mol objects.
"""

import numpy as np
import scipy.cluster.hierarchy as sch
from rdkit.Chem import Mol
from rdkit.Chem.rdFingerprintGenerator import FingerprintGenerator64, GetMorganGenerator
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
from scipy.spatial.distance import euclidean


def _get_fingerprints(mol_list: list[Mol]) -> list[FingerprintGenerator64]:
    """Generate Morgan fingerprints for a list of RDKit molecules.

    Args:
        mol_list (list[Mol]): List of RDKit molecules.

    Returns:
        list[FingerprintGenerator64]: List of Morgan fingerprints.
    """
    return [
        GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(mol)
        for mol in mol_list
    ]


def hierarchical_clustering(mol_list: list, cluster_method: str, t: None | int | float):
    """Perform hierarchical clustering on a list of RDKit molecules using Tanimoto similarity, with logging."""
    yield {"log": "🔍 Generating fingerprints..."}
    fingerprints = _get_fingerprints(mol_list)

    yield {"log": "🧮 Computing pairwise Tanimoto similarity matrix..."}
    n_mols = len(mol_list)
    sim_matrix = np.zeros((n_mols, n_mols))

    for i in range(n_mols):
        for j in range(i + 1, n_mols):
            sim = TanimotoSimilarity(fingerprints[i], fingerprints[j])
            sim_matrix[i, j] = sim
            sim_matrix[j, i] = sim

    yield {"log": "📏 Converting similarity to distance matrix..."}
    dist_matrix = 1 - sim_matrix

    yield {"log": "🌿 Performing hierarchical clustering..."}
    linkage_matrix = sch.linkage(dist_matrix, method="average")

    yield {"log": "🔗 Creating cluster labels..."}
    cluster_labels = sch.fcluster(linkage_matrix, criterion=cluster_method, t=t)

    yield {"result": cluster_labels}


def find_cluster_centroids(mol_list: list, cluster_labels):
    """Find centroids of clusters based on Euclidean distance from mean fingerprint, with logging."""
    yield {"log": "📍 Calculating fingerprints for centroid selection..."}
    fingerprints = _get_fingerprints(mol_list)

    unique_clusters = set(cluster_labels)
    centroids = {}

    for cluster_id in unique_clusters:
        yield {"log": f"  🔸 Processing cluster {cluster_id}..."}

        cluster_indices = [
            i for i, label in enumerate(cluster_labels) if label == cluster_id
        ]
        cluster_fingerprints = [fingerprints[i] for i in cluster_indices]

        # Convert to NumPy array
        fp_matrix = np.array(cluster_fingerprints)
        centroid_vector = np.mean(fp_matrix, axis=0)

        distances = [euclidean(fp, centroid_vector) for fp in fp_matrix]
        closest_idx = np.argmin(distances)
        centroids[cluster_id] = mol_list[cluster_indices[closest_idx]]

        yield {
            "log": f"    ✅ Cluster {cluster_id} centroid selected (index {cluster_indices[closest_idx]})."
        }

    yield {"result": centroids}
