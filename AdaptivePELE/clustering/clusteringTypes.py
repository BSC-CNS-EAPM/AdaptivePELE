from AdaptivePELE.constants import blockNames


class CLUSTERING_TYPES:
    contacts, contactMapAccumulative, lastSnapshot = range(3)

CLUSTERING_TYPE_TO_STRING_DICTIONARY = {
    CLUSTERING_TYPES.contacts: blockNames.ClusteringTypes.contacts,
    CLUSTERING_TYPES.contactMapAccumulative: blockNames.ClusteringTypes.contactMapAccumulative,
    CLUSTERING_TYPES.lastSnapshot: blockNames.ClusteringTypes.lastSnapshot
}


class SIMILARITY_TYPES:
    differenceDistance, Jaccard, correlation = range(3)

SIMILARITY_TYPES_TO_STRING_DICTIONARY = {
    SIMILARITY_TYPES.differenceDistance: blockNames.ClusteringTypes.differenceDistance,
    SIMILARITY_TYPES.Jaccard: blockNames.ClusteringTypes.Jaccard,
    SIMILARITY_TYPES.correlation: blockNames.ClusteringTypes.correlation
}
