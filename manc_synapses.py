import neuprint
import os
import pandas as pd
client = neuprint.Client('neuprint-pre.janelia.org', dataset='vnc')

def fetch_downstream_connections(neuron_id):
    t1_ROIs = ["IntTct","LTct","LegNp(T1)(L)","LegNp(T1)(R)","NTct(UTct-T1)(L)","NTct(UTct-T1)(R)","mVAC(T1)(L)","mVAC(T1)(R)"]
    t1_area = neuprint.queries.NeuronCriteria(rois=["LegNp(T1)(L)","LegNp(T1)(R)"], roi_req="any")
    connections = neuprint.fetch_simple_connections([neuron_id], downstream_criteria=t1_area, rois=t1_ROIs)
    synapse_values = connections[["bodyId_post", "weight"]].copy()
    synapse_values = synapse_values.set_index("bodyId_post")
    synapse_values = synapse_values.loc[synapse_values['weight'] > 10]
    # find neurons that go into the abdomen and remove them
    ROIs_to_avoid = ["ANm","LegNp(T3)(L)","LegNp(T3)(R)","HTct(UTct-T3)(L)","HTct(UTct-T3)(R)"]
    criteria = neuprint.queries.NeuronCriteria(bodyId=synapse_values.index.to_list(), rois=ROIs_to_avoid, roi_req="any")
    neurons_to_remove,_ = neuprint.fetch_neurons(criteria)
    neurons_to_remove = neurons_to_remove.bodyId.to_list()
    return synapse_values.loc[~synapse_values.index.isin(neurons_to_remove)]

def get_percent_input(neuron_id):
    conn_table = fetch_downstream_connections(neuron_id)
    IDs = conn_table.index.to_list()
    if not len(IDs):
        return pd.DataFrame({"percent":[], "weight":[]})
    neurons,_ = neuprint.fetch_neurons(IDs)
    neurons = neurons.set_index("bodyId")
    neurons = neurons.reindex(IDs)
    # remove fragments by filtering out IDs with no soma location
    dendrite_number = neurons[~neurons["somaLocation"].isnull()]
    dendrite_number = dendrite_number["post"]
    inputs = conn_table.merge(dendrite_number, left_index=True, right_index=True, how="right")
    inputs["percent"] = (inputs["weight"]/inputs["post"])*100
    return inputs.sort_values("percent",ascending=False)

def fetch_upstream_connections(neuron_id):
    connections = neuprint.fetch_simple_connections(None,[neuron_id])
    synapse_values = connections[["bodyId_pre", "weight"]].copy()
    synapse_values = synapse_values.set_index("bodyId_pre")
    return synapse_values

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

def cascade_csvs(start_neuron, threshold):
    def make_csvs_in_list(neuron_list):
        downstream_neurons = []
        for neuron_id in neuron_list:
            print(neuron_id, end=" ")
            if str(neuron_id)+"_downstreampartners.csv" in files:
                print("file already exists")
                pass
            else:
                downstream_new = make_csv(neuron_id, threshold, folder)
                downstream_neurons = downstream_neurons+downstream_new
                if len(downstream_new):
                    print("created file")
                else:
                    print("no downstream partners")
        return downstream_neurons
    # takes an initial starting neuron and finds its most significant downstream partners
    # then does the same for each of the downstream partners for n layers
    folder = str(start_neuron)+"-"+str(threshold)+"/"
    if not os.path.exists(folder):
        os.makedirs(folder)
    files = os.listdir(folder)
    print("downstream of",start_neuron)
    first_layer = make_csvs_in_list([start_neuron])
    print("second layer")
    threshold=threshold
    second_layer = make_csvs_in_list(first_layer)
    print("third_layer")
    threshold=threshold
    third_layer = make_csvs_in_list(second_layer)
    #print("fourth layer")
    #make_csvs_in_list(third_layer)
    
def make_csv(neuron_id, threshold, folder=""):
    dataframe = get_percent_input(neuron_id)
    dataframe = dataframe.loc[dataframe['percent'] >= threshold]
    #dataframe = dataframe.loc[dataframe['weight'] >= threshold]
    dataframe = dataframe.sort_values("weight",ascending=False)
    if dataframe.size == 0:
        return []
    dataframe.to_csv(folder+str(neuron_id)+"_downstreampartners.csv")
    return dataframe.index.to_list()

def synapses_between(upstream,downstream):
    return neuprint.fetch_simple_connections([upstream], [downstream])