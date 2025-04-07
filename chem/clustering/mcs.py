import numpy as np
from rdkit.Chem import rdFMCS
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform


def compute_mcs_similarity(mol1, mol2):
    """Compute similarity based on Maximum Common Substructure (MCS)."""
    mcs = rdFMCS.FindMCS([mol1, mol2])
    return mcs.numAtoms  # Larger MCS = More similar structures


def pairwise_mcs_distance(fragments):
    """Compute pairwise MCS-based distance (1 - normalized MCS size)."""
    n = len(fragments)
    distance_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            mcs_size = compute_mcs_similarity(fragments[i], fragments[j])
            max_atoms = max(fragments[i].GetNumAtoms(), fragments[j].GetNumAtoms())
            distance = 1 - (mcs_size / max_atoms)  # Normalize distance
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance  # Symmetric matrix

    return distance_matrix


def hierarchical_clustering(mol_list, cluster_method, t):
    """Cluster fragments using hierarchical clustering based on MCS similarity."""
    dist_matrix = pairwise_mcs_distance(mol_list)
    condensed_matrix = squareform(
        dist_matrix
    )  # Convert to condensed form for clustering

    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_matrix, method="ward")

    # Associate fragments to clusters. Indexes are the fgragments positions, value is the cluster number
    cluster_labels = fcluster(linkage_matrix, criterion=cluster_method, t=t)

    return cluster_labels


def find_cluster_centroids(mol_list, cluster_labels):
    """Find the centroid fragment for each cluster (smallest average distance to others)."""
    unique_clusters = set(cluster_labels)
    centroids = {}

    for cluster_id in unique_clusters:
        cluster_indices = [
            i for i, label in enumerate(cluster_labels) if label == cluster_id
        ]
        sub_matrix = pairwise_mcs_distance([mol_list[i] for i in cluster_indices])

        # Find the fragment with the smallest average distance
        avg_distances = np.mean(sub_matrix, axis=1)
        centroid_idx = cluster_indices[np.argmin(avg_distances)]
        centroids[cluster_id] = mol_list[centroid_idx]

    return centroids
