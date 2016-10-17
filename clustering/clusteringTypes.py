from constants import blockNames


class CLUSTERING_TYPES:
    contacts, contactMapAffinity, contactMapAgglomerative, contactMapAccumulative = range(4)

CLUSTERING_TYPE_TO_STRING_DICTIONARY = {
    CLUSTERING_TYPES.contacts: blockNames.ClusteringTypes.contacts,
    CLUSTERING_TYPES.contactMapAffinity: blockNames.ClusteringTypes.contactMapAffinity,
    CLUSTERING_TYPES.contactMapAgglomerative: blockNames.ClusteringTypes.contactMapAgglomerative,
    CLUSTERING_TYPES.contactMapAccumulative: blockNames.ClusteringTypes.contactMapAccumulative
}


class SIMILARITY_TYPES:
    differenceDistance, Jaccard, correlation = range(3)

SIMILARITY_TYPES_TO_STRING_DICTIONARY = {
    SIMILARITY_TYPES.differenceDistance: blockNames.ClusteringTypes.differenceDistance,
    SIMILARITY_TYPES.Jaccard: blockNames.ClusteringTypes.Jaccard,
    SIMILARITY_TYPES.correlation: blockNames.ClusteringTypes.correlation
}
