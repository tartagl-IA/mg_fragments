"""Clustering module for cheminformatics."""

from enum import Enum


class ClusteringMethod(Enum):
    """Enum for clustering methods."""

    DIST = "distance"
    MAX_CLUSTERS = "maxclust"


class ClusteringType(Enum):
    """Enum for clustering types."""

    MCS = "Maximum Common Substructure (MCS)"
    TANIMOTO = "Tanimoto Similarity"


__all__ = [
    "ClusteringMethod",
    "ClusteringType",
]
