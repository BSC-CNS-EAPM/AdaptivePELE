import blockNames

class CLUSTERING_TYPES:
    contacts, contactMapAffinity, contactMapAgglomerative = range(3)

CLUSTERING_TYPE_TO_STRING_DICTIONARY = {
    CLUSTERING_TYPES.contacts : blockNames.ClusteringTypes.contacts,
    CLUSTERING_TYPES.contactMapAffinity : blockNames.ClusteringTypes.contactMapAffinity,
    CLUSTERING_TYPES.contactMapAgglomerative : blockNames.ClusteringTypes.contactMapAgglomerative
}

