import neuprint
import os

def fetch_downstream_synapses(neuron_id):
    client = neuprint.Client('neuprint.janelia.org', dataset='manc:v1.0')
    synapses = neuprint.fetch_synapse_connections(source_criteria=neuron_id, batch_size=100)
    syn = synapses[["bodyId_post", "x_pre", "y_pre", "z_pre"]].copy()
    return syn

def fetch_upstream_synapses(neuron_id):
    client = neuprint.Client('neuprint.janelia.org', dataset='manc:v1.0')
    synapses = neuprint.fetch_synapse_connections(target_criteria=neuron_id)
    syn = synapses[["bodyId_pre", "x_pre", "y_pre", "z_pre"]].copy()
    return syn

def synapse_y_limit(synapses, y_limit):
    return synapses[synapses["y_pre"]<y_limit]

def get_value_counts(df):
    return df.iloc[:,0].value_counts()

def downstream_of(neuron_id, N):
    #get synapses and sort value counts up to the Nth downstream synapse
    downstream_synapses = fetch_downstream_synapses(neuron_id)
    values = get_value_counts(downstream_synapses)
    indices = values.head(N).index.to_list()
    for index in indices:
        print(index)
    synapse_values = values.head(N).to_list()
    for synapse_value in synapse_values:
        print(synapse_value)

def get_type(neuron_ids):
    client = neuprint.Client('neuprint.janelia.org', dataset='manc:v1.0')
    neurons, _ = neuprint.fetch_neurons(neuron_ids)
    neurons = neurons[["bodyId", "type"]].copy()
    neurons = neurons.set_index("bodyId")
    neurons = neurons.reindex(neuron_ids)
    return neurons

def print_type(N):
    IDs = []
    for _ in range(N):
        ID = input()
        IDs.append(int(ID))
    types = get_type(IDs)
    print(types)
    for t in types["type"].items():
        print(t[1])

def cascade_csvs(neuron_id, N):
    # takes an initial starting neuron and finds its most significant downstream partners
    # then does the same for each of the downstream partners for 3 layers
    files = os.listdir()
    downstream_neurons = make_csv(neuron_id, N)
    for neuron_id in downstream_neurons:
        print(neuron_id, end=" ")
        if str(neuron_id)+"_downstreampartners.csv" in files:
            print("file already exists")
            pass
        else:
            make_csv(neuron_id, N)
            print("created file")
    


def make_csv(neuron_id, N):
    downstream_synapses = fetch_downstream_synapses(neuron_id)
    values = get_value_counts(downstream_synapses).head(N)
    types = get_type(values.index.to_list())
    dataframe = types.merge(values, left_index=True, right_index=True)
    dataframe.to_csv(str(neuron_id)+"_downstreampartners.csv")
    return values.index.to_list()