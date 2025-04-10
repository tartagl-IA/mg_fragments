"""Minimum Common Substructure (MCS) clustering for molecular fragments.

This module provides functions to compute MCS-based similarity
and perform hierarchical clustering on molecular fragments.
It also includes functions to find cluster centroids based on MCS similarity.
The MCS similarity is computed using the RDKit library.
The clustering is performed using SciPy's hierarchical clustering methods.
The module is designed to work with molecular structures represented as RDKit Mol objects.
"""

import time
from typing import Generator

import numpy as np
from rdkit.Chem import Mol, rdFMCS
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform


def compute_mcs_similarity(mol1: Mol, mol2: Mol) -> int:
    """Compute MCS similarity between two molecules.

    Args:
        mol1 (Mol): First molecule.
        mol2 (Mol): Second molecule.

    Returns:
        int: Number of atoms in the MCS.
    """
    mcs = rdFMCS.FindMCS([mol1, mol2])
    return mcs.numAtoms  # Larger MCS = More similar structures


def compute_mcs_similarity_timeout(mol1: Mol, mol2: Mol, timeout: int = 5) -> int:
    """Compute MCS similarity with a timeout.

    Args:
        mol1 (Mol): First molecule.
        mol2 (Mol): Second molecule.
        timeout (int): Timeout in seconds.

    Returns:
        int: Number of atoms in the MCS or 0 if timeout occurs.
    """
    try:
        res = rdFMCS.FindMCS([mol1, mol2], timeout=timeout)
        return res.numAtoms  # Larger MCS = More similar structures
    except Exception:
        return 0


def pairwise_mcs_distance(fragments: list[Mol]) -> np.ndarray:
    """Compute pairwise MCS-based distance (1 - normalized MCS size)."""
    n = len(fragments)
    distance_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            mcs_size = compute_mcs_similarity_timeout(fragments[i], fragments[j])
            max_atoms = max(fragments[i].GetNumAtoms(), fragments[j].GetNumAtoms())
            distance = 1 - (mcs_size / max_atoms)  # Normalize distance
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance  # Symmetric matrix

    return distance_matrix


def hierarchical_clustering(
    mol_list: list[Mol], cluster_method: str, t: None | int | float
) -> Generator[dict[str, str | dict], None, None]:
    """Cluster fragments using hierarchical clustering based on MCS similarity."""
    yield {"log": "🔍 Starting pairwise MCS distance computation..."}
    start = time.time()
    dist_matrix = pairwise_mcs_distance(mol_list)
    yield {"log": f"✅ pairwise_mcs_distance done in {time.time() - start:.2f}s"}

    yield {"log": "📐 Converting to condensed matrix..."}
    condensed_matrix = squareform(dist_matrix)

    yield {"log": "🌿 Performing hierarchical clustering..."}
    linkage_matrix = linkage(condensed_matrix, method="ward")

    yield {"log": "🔗 Assigning cluster labels..."}
    cluster_labels = fcluster(linkage_matrix, criterion=cluster_method, t=t)

    # Final result: yield with a special key
    yield {"result": cluster_labels}


def find_cluster_centroids(
    mol_list, cluster_labels
) -> Generator[dict[str, str | dict], None, None]:
    """Find the centroid fragment for each cluster (smallest average distance to others), yielding logs."""
    yield {"log": "📍 Finding cluster centroids..."}
    unique_clusters = set(cluster_labels)
    centroids = {}

    for cluster_id in unique_clusters:
        yield {"log": f"🔸 Processing cluster {cluster_id}..."}

        cluster_indices = [
            i for i, label in enumerate(cluster_labels) if label == cluster_id
        ]
        cluster_mols = [mol_list[i] for i in cluster_indices]

        start = time.time()
        sub_matrix = pairwise_mcs_distance(cluster_mols)
        duration = time.time() - start
        yield {
            "log": f"⏱️ MCS distance matrix computed in {duration:.2f}s for {len(cluster_indices)} fragments."
        }

        avg_distances = np.mean(sub_matrix, axis=1)
        centroid_idx = cluster_indices[np.argmin(avg_distances)]
        centroids[cluster_id] = mol_list[centroid_idx]

        yield {
            "log": f"✅ Cluster {cluster_id} centroid selected (index {centroid_idx})."
        }

    yield {"result": centroids}
