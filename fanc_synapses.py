import caveclient as cv
import datetime

def fetch_downstream_synapses(neuron_id):
    client = cv.CAVEclient('fanc_production_mar2021')
    synapse_table = client.info.get_datastack_info()['synapse_table']
    synapses = client.materialize.live_query(synapse_table,datetime.datetime.now(datetime.timezone.utc),filter_equal_dict = {'pre_pt_root_id': neuron_id})
    syn = synapses[["post_pt_root_id", "pre_pt_position"]].copy()
    return syn

def fetch_upstream_synapses(neuron_id):
    client = cv.CAVEclient('fanc_production_mar2021')
    synapse_table = client.info.get_datastack_info()['synapse_table']
    synapses = client.materialize.live_query(synapse_table,datetime.datetime.now(datetime.timezone.utc),filter_equal_dict = {'post_pt_root_id': neuron_id})
    syn = synapses[["pre_pt_root_id", "pre_pt_position"]].copy()
    return syn

def synapse_y_limit(synapses, y_limit):
    y_axis = []
    for _, value in synapses["pre_pt_position"].items():
        y_axis.append(value[1])
    synapses["pre_pt_y_axis"] = y_axis
    return synapses[synapses["pre_pt_y_axis"]<y_limit]

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

