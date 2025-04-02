import numpy as np
import scipy.cluster.hierarchy as sch
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
from scipy.spatial.distance import euclidean


def _get_fingerprints(mol_list):
    return [
        GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(mol)
        for mol in mol_list
    ]


def hierarchical_clustering(mol_list, cluster_method, t):
    fingerprints = _get_fingerprints(mol_list)

    # Compute pairwise similarity matrix
    n_mols = len(mol_list)
    sim_matrix = np.zeros((n_mols, n_mols))

    for i in range(n_mols):
        for j in range(i + 1, n_mols):
            sim_matrix[i, j] = TanimotoSimilarity(fingerprints[i], fingerprints[j])
            sim_matrix[j, i] = sim_matrix[i, j]  # Symmetric matrix

    # Convert similarity to distance (1 - similarity)
    dist_matrix = 1 - sim_matrix

    # Perform hierarchical clustering
    linkage_matrix = sch.linkage(dist_matrix, method="average")

    # Create clusters
    cluster_labels = sch.fcluster(linkage_matrix, criterion=cluster_method, t=t)

    return cluster_labels


def find_cluster_centroids(mol_list, cluster_labels):

    fingerprints = _get_fingerprints(mol_list)

    unique_clusters = set(cluster_labels)
    centroids = {}

    for cluster_id in unique_clusters:
        cluster_indices = [
            i for i, label in enumerate(cluster_labels) if label == cluster_id
        ]
        cluster_fingerprints = [fingerprints[i] for i in cluster_indices]
        centroid = np.mean(cluster_fingerprints, axis=0)
        distances = [euclidean(fp, centroid) for fp in cluster_fingerprints]
        closest_idx = np.argmin(distances)
        centroids[cluster_id] = mol_list[cluster_indices[closest_idx]]

    return centroids
