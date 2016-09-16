import blockNames

class CLUSTERING_TYPES:
    contacts, contactMap, agglomerative = range(3)

CLUSTERING_TYPE_TO_STRING_DICTIONARY = {
    CLUSTERING_TYPES.contacts : blockNames.ClusteringTypes.contacts,
    CLUSTERING_TYPES.contactMap : blockNames.ClusteringTypes.contactMap,
    CLUSTERING_TYPES.agglomerative : blockNames.ClusteringTypes.agglomerative
}

