import blockNames

class CLUSTERING_TYPES:
    contacts, contactMap = range(2)

CLUSTERING_TYPE_TO_STRING_DICTIONARY = {
    CLUSTERING_TYPES.contacts : blockNames.ClusteringTypes.contacts,
    CLUSTERING_TYPES.contactMap : blockNames.ClusteringTypes.contactMap
}

class THRESHOLD_CALCULATOR_TYPES:
    heaviside, constant = range(2)

THRESHOLD_CALCULATOR_TYPE_TO_STRING_DICTIONARY = {
    THRESHOLD_CALCULATOR_TYPES.heaviside : blockNames.ThresholdCalculator.heaviside,
    THRESHOLD_CALCULATOR_TYPES.constant : blockNames.ThresholdCalculator.constant
}
