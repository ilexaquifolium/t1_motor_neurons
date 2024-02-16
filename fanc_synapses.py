import caveclient as cv
import datetime
import pandas as pd
import os
import seaserpent as ss
client = cv.CAVEclient('fanc_production_mar2021')

def fetch_downstream_synapses(neuron_id):
    synapses = client.materialize.synapse_query(pre_ids=[neuron_id])
    syn = synapses[["post_pt_root_id", "pre_pt_position"]].copy()
    return syn

def get_percent_input(conn_table):
    IDs = conn_table.index.to_list()
    synapses = client.materialize.synapse_query(post_ids=IDs)
    dendrite_number = synapses["post_pt_root_id"].value_counts().rename("post")
    inputs = pd.concat([conn_table,dendrite_number], axis=1)
    inputs["percent"] = (inputs["count"]/inputs["post"])*100
    return inputs.sort_values("percent",ascending=False)

def fetch_upstream_synapses(neuron_id):
    synapses = client.materialize.synapse_query(post_ids=[neuron_id])
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

def downstream_of(neuron_id, threshold):
    #neuron_id = yesterday(neuron_id)
    #get downstream neurons and sort by synapse counts up to the Nth downstream synapse
    downstream_synapses = fetch_downstream_synapses(neuron_id)
    t1_downstream_synapses = synapse_y_limit(downstream_synapses,118000)
    values = get_value_counts(t1_downstream_synapses)
    values = values.loc[values >= threshold]
    percentage_table = get_percent_input(values)
    # filter out fragments
    percentage_table = percentage_table.loc[percentage_table.post > 100]
    return percentage_table

def get_type(neuron_ids):
    neurons = pd.read_csv("fanc_types.csv", index_col=0)
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
    # then does the same for each of the downstream partners for 3 layers
    folder = str(start_neuron)+"-"+str(threshold)+"/"
    if not os.path.exists(folder):
        os.makedirs(folder)
    files = os.listdir(folder)
    print("downstream of",start_neuron)
    first_layer = make_csvs_in_list([start_neuron])
    print("second layer")
    #threshold=1*threshold
    second_layer = make_csvs_in_list(first_layer)
    print("third_layer")
    #threshold=1*threshold
    third_layer = make_csvs_in_list(second_layer)
    #print("fourth layer")
    #make_csvs_in_list(third_layer)
    
def make_csv(neuron_id, percentage_threshold=0.5, folder=""):
    table = downstream_of(neuron_id,3)
    dataframe = table.loc[table.percent >= percentage_threshold]
    if dataframe.size == 0:
        return []
    dataframe = dataframe.rename(columns={'count': 'weight'})
    dataframe.index.names = ['bodyId']
    dataframe.to_csv(folder+str(neuron_id)+"_downstreampartners.csv")
    return dataframe.index.to_list()

def update_fanc_id(fanc_id):
    if fanc_id == "" or fanc_id == None or fanc_id == "NotAssigned":
        return None
    # this produces a FutureWarning
    newest_id = client.chunkedgraph.suggest_latest_roots(fanc_id)
    return newest_id

def mf_match(manc_id, seatable=None):
    if not seatable:
        seatable = ss.Table(table='fanc851_manc_nblast95_60')
    row = seatable[seatable.queryID == str(manc_id)].iloc[0]
    if row.manualAssignment:
        return update_fanc_id(row.manualAssignment)
    try:
        nblast_match = row.nBlastMatchID
    except:
        nblast_match = None
    try:
        cosine_match = row.conMatchID
    except:
        cosine_match = None
    if nblast_match == cosine_match:
        return update_fanc_id(nblast_match)
    else:
        if row.check_L_R[0] == "Yes":
            symmetry_string = ", also neuron may be left/right flipped"
        else:
            symmetry_string = ""
        print("nblast and cosine similarity don't agree on a match for ", manc_id, symmetry_string)
        return ("nblast:", update_fanc_id(nblast_match), "connectivity:", update_fanc_id(cosine_match))

def at_time_of_match(fanc_id):
    return client.chunkedgraph.suggest_latest_roots(fanc_id, timestamp=datetime.datetime(2023,7,20,8,10,1))

def fm_match(fanc_id, seatable=None):
    if not seatable:
        seatable = ss.Table(table='fanc851_manc_nblast95_60')
    # search in the manual matches
    past_ids = [str(x) for x in client.chunkedgraph.get_past_ids(fanc_id)["past_id_map"][fanc_id]]
    past_ids.append(str(fanc_id))
    past_ids.append(str(at_time_of_match(fanc_id)))
    rows = seatable[seatable.manualAssignment.isin(past_ids)]
    if rows.shape[0] == 1:
        print("manually assigned")
        return int(rows.iloc[0].queryID)

    # if there are no manual matches
    if rows.shape[0] == 0:
        # get equivalent ID from time when matching was done
        fanc_id = at_time_of_match(fanc_id)
        rows = seatable[seatable.nBlastMatchID.isin(past_ids) or seatable.conMatchID.isin(past_ids)]
    
    # if there's one row and the nblast and connection matches are the same
    if rows.shape[0] == 1 and rows.iloc[0].nBlastMatchID == rows.iloc[0].conMatchID:
        return int(rows.iloc[0].queryID)
    else:
        print("no definite match for", fanc_id, ", options returned as a list.")
    return rows.queryID.to_list()

def rank_matches(fanc_csv, manc_csv):
    table = ss.Table(table='fanc851_manc_nblast95_60')
    fanc = pd.read_csv(fanc_csv, index_col=0)
    manc = pd.read_csv(manc_csv, index_col=0)
    #fanc = fanc_csv
    #manc = manc_csv
    assigned_fanc = []
    clean_fanc = []
    for id, _ in manc.iterrows():
        fanc_id = mf_match(id,table)
        assigned_fanc.append(fanc_id)
        if not isinstance(fanc_id, tuple):
            clean_fanc.append(fanc_id)
        else:
            clean_fanc.append(None)
    manc["fanc"] = clean_fanc
    fanc.index = fanc.index.astype(float)
    #return manc, fanc
    combined_table = manc.merge(fanc, left_on="fanc", right_index=True, how="outer")
    #combined_table["fanc"] = assigned_fanc
    return combined_table

def rank_f_f(fanc1, fanc2):
    fanc1 = pd.read_csv(fanc1, index_col=0)
    fanc2 = pd.read_csv(fanc2, index_col=0)
    combined_table = fanc1.merge(fanc2, left_index=True, right_index=True, how="outer")
    return combined_table

def yesterday(id):
    return client.chunkedgraph.suggest_latest_roots(id, timestamp=(datetime.datetime.now() - datetime.timedelta(days = 1)))

def synapses_between(upstream,downstream):
    upstream = yesterday(upstream)
    downstream = yesterday(downstream)
    return client.materialize.synapse_query(pre_ids=[upstream],post_ids=[downstream])

def update_MN_sheet():
    # update the column of current FANC IDs on the T1 Motor Neuron matching sheet
    csv_filename = "FANC-MANC-matching - MNs.csv"
    T1MNs = pd.read_csv(csv_filename)
    latest_fanc_id = []
    for n, fanc_id in enumerate(T1MNs.fanc_root_id):
        print(n,"/",T1MNs.shape[0],end=", ", flush=True)
        latest_fanc_id.append(client.chunkedgraph.suggest_latest_roots(fanc_id))
    T1MNs["latest_fanc_id"] = latest_fanc_id
    T1MNs.to_csv(csv_filename)
    print("saved to",csv_filename)

def MN_match(neuron_id):
    csv_filename = "FANC-MANC-matching - MNs.csv"
    T1MNs = pd.read_csv(csv_filename)
    if len(str(neuron_id)) < 12:#manc id
        return T1MNs.loc[T1MNs.manc_bodyid == neuron_id].latest_fanc_id
    else:#fanc id
        neuron_id = update_fanc_id(neuron_id)
        return T1MNs.loc[T1MNs.latest_fanc_id == neuron_id].manc_bodyid
