import neuprint
import os

def fetch_downstream_connections(neuron_id):
    client = neuprint.Client('neuprint.janelia.org', dataset='manc:v1.0')
    connections = neuprint.fetch_simple_connections([neuron_id])
    synapse_values = connections[["bodyId_post", "weight"]].copy()
    synapse_values = synapse_values.set_index("bodyId_post")
    return synapse_values

def fetch_upstream_connections(neuron_id):
    client = neuprint.Client('neuprint.janelia.org', dataset='manc:v1.0')
    connections = neuprint.fetch_simple_connections(None,[neuron_id])
    synapse_values = connections[["bodyId_pre", "weight"]].copy()
    synapse_values = synapse_values.set_index("bodyId_pre")
    return synapse_values

def get_value_counts(df):
    return df.iloc[:,0].value_counts()

def downstream_of(neuron_id, threshold):
    #get downstream neurons and sort by synapse counts up to the Nth downstream synapse
    values = fetch_downstream_connections(neuron_id)
    values = values.loc[values['weight'] >= threshold]
    indices = values.index.to_list()
    for index in indices:
        print(index)
    synapse_values = values.weight.to_list()
    for synapse_value in synapse_values:
        print(synapse_value)

def get_type(neuron_ids):
    client = neuprint.Client('neuprint.janelia.org', dataset='manc:v1.0')
    neurons, _ = neuprint.fetch_neurons(neuron_ids)
    neurons = neurons[["bodyId", "type"]].copy()
    neurons = neurons.set_index("bodyId")
    neurons = neurons.reindex(neuron_ids)
    return neurons

def print_type():
    IDs = []
    ID = input()
    while ID != "":
        IDs.append(int(ID))
        ID = input()
    types = get_type(IDs)
    print(types)
    for t in types["type"].items():
        print(t[1])

def cascade_csvs(neuron_id, threshold):
    # takes an initial starting neuron and finds its most significant downstream partners
    # then does the same for each of the downstream partners for 3 layers
    files = os.listdir()
    downstream_neurons = make_csv(neuron_id, threshold)
    for neuron_id in downstream_neurons:
        print(neuron_id, end=" ")
        if str(neuron_id)+"-downstreampartners.csv" in files:
            print("file already exists")
            pass
        else:
            make_csv(neuron_id, threshold)
            print("created file")
    
def make_csv(neuron_id, threshold):
    values = fetch_downstream_connections(neuron_id)
    values = values.loc[values['weight'] >= threshold]
    if values.size == 0:
        return []
    types = get_type(values.index.to_list())
    dataframe = types.merge(values, left_index=True, right_index=True)
    dataframe.to_csv(str(neuron_id)+"_downstreampartners.csv")
    return values.index.to_list()