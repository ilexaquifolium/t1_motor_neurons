import caveclient as cv
import datetime
import pandas as pd
import os
import seaserpent as ss

def fetch_downstream_synapses(neuron_id, client=None):
    """Return a dataframe of synapses that are downstream of the neuron.

    Parameters
    ----------
    neuron_id :         int or str
                        FANC neuron ID to find the downstream synapses of
    client :            caveclient.frameworkclient.CAVEclientFull
                        CAVEclient to query information from. If not passed to the function it will create its own.

    Returns
    -------
    pandas DataFrame
                        DataFrame with columns -
                            post_pt_root_id: the ID of the downstream neuron
                            pre_pt_position: the [x,y,z] coordinates in list format of the presynapse

    """
    if not client:
        client = cv.CAVEclient('fanc_production_mar2021')
    synapses = client.materialize.synapse_query(pre_ids=[neuron_id])
    syn = synapses[["post_pt_root_id", "pre_pt_position"]].copy()
    return syn

def get_percent_input(conn_table, client=None):
    """Take a Series of synapse counts and return a dataframe of neurons with net and percentage inputs.

    Parameters
    ----------
    conn_table :        pandas Series
                        index is the ID of the downstream neuron, values are the number of connections to this neuron
    client :            caveclient.frameworkclient.CAVEclientFull
                        CAVEclient to query information from. If not passed to the function it will create its own.

    Returns
    -------
    pandas DataFrame
                        DataFrame sorted by high to low percent with columns -
                            post_pt_root_id: the ID of the downstream neuron
                            count: the number of connections to the downstream neuron
                            inputs: the total number of inputs to the downstream neuron
                            percent: the percentage of input to the downstream
                                neuron represented by the connections in the input

    """
    if not client:
        client = cv.CAVEclient('fanc_production_mar2021')
    IDs = conn_table.index.to_list()
    synapses = client.materialize.synapse_query(post_ids=IDs)
    dendrite_number = synapses["post_pt_root_id"].value_counts().rename("inputs")
    inputs = pd.concat([conn_table,dendrite_number], axis=1)
    inputs["percent"] = (inputs["count"]/inputs["inputs"])*100
    return inputs.sort_values("percent",ascending=False)

def fetch_upstream_synapses(neuron_id, client=None):
    """Return a dataframe of synapses that are upstream of the neuron.

    Parameters
    ----------
    neuron_id :         int or str
                        FANC neuron ID to find the upstream synapses of
    client :            caveclient.frameworkclient.CAVEclientFull
                        CAVEclient to query information from. If not passed to the function it will create its own.

    Returns
    -------
    pandas DataFrame
                        DataFrame with columns -
                            pre_pt_root_id: the ID of the upstream neuron
                            pre_pt_position: the [x,y,z] coordinates in list format of the presynapse

    """
    if not client:
        client = cv.CAVEclient('fanc_production_mar2021')
    synapses = client.materialize.synapse_query(post_ids=[neuron_id])
    syn = synapses[["pre_pt_root_id", "pre_pt_position"]].copy()
    return syn

def synapse_y_limit(synapses, posterior_limit, anterior_limit=0):
    """Take a dataframe of synapses and return only those which are above the posterior threshold
    and below the anterior threshold in the VNC.

    Parameters
    ----------
    pandas DataFrame
                        DataFrame with columns -
                            post_pt_root_id: the ID of the downstream neuron
                            pre_pt_position: the [x,y,z] coordinates in list format of the presynapse
    posterior_limit :   int
                        y value in the VNC below which synapses will be discarded
    anterior_limit :   int
                        y value in the VNC above which synapses will be discarded

    Returns
    -------
    pandas DataFrame
                        DataFrame with columns -
                            post_pt_root_id/pre_pt_root_id: the ID of the downstream neuron
                            pre_pt_position: the [x,y,z] coordinates of the presynapse in list format
                            pre_pt_y_axis: just the y value of the coordinates

    """
    y_axis = []
    for _, value in synapses["pre_pt_position"].items():
        y_axis.append(value[1])
    synapses["pre_pt_y_axis"] = y_axis
    return synapses[synapses["pre_pt_y_axis"].between(anterior_limit,posterior_limit)]

def downstream_of(neuron_id, threshold=3, client=None):
    """Return a dataframe of neurons that are downstream of the neuron, by percentage input.

    Parameters
    ----------
    neuron_id :         int or str
                        FANC neuron ID to find the downstream partners of
    threshold :         int
                        minimum number of synapses required to be considered a connection between neurons
    client :            caveclient.frameworkclient.CAVEclientFull
                        CAVEclient to query information from. If not passed to the function it will create its own.

    Returns
    -------
    pandas DataFrame
                        DataFrame sorted by high to low percent with columns -
                            post_pt_root_id: the ID of the downstream neuron
                            count: the number of connections to the downstream neuron
                            inputs: the total number of inputs to the downstream neuron
                            percent: the percentage of input to the downstream
                                neuron represented by the connections in the input

    """
    if not client:
        client = cv.CAVEclient('fanc_production_mar2021')
    #get downstream neurons and sort by synapse counts up to the Nth downstream synapse
    downstream_synapses = fetch_downstream_synapses(neuron_id, client)
    t1_downstream_synapses = synapse_y_limit(downstream_synapses,118000)

    # get value counts as Series
    values = t1_downstream_synapses.iloc[:,0].value_counts()
    values = values.loc[values >= threshold]
    percentage_table = get_percent_input(values, client)

    # filter out fragments
    percentage_table = percentage_table.loc[percentage_table.inputs > 100]
    return percentage_table

def cascade_csvs(start_neuron, percentage_threshold=1, connection_threshold=3, layers=3, client=None):
    """Take an initial starting neuron and save .csv files of its most significant downstream partners
        then do the same for each of the downstream partners for the chosen number of layers. The dataframes
        representing the downstream partners of each neuron are saved in a standard .csv
        format inside a folder named after the start neuron

    Parameters
    ----------
    start_neuron :          int
                            FANC neuron ID to find the downstream partners of
    percentage_threshold :  float
                            minimum percentage input required to be included in the saved csv files
    connection_threshold :  int
                            minimum number of synapses required to be considered a connection between neurons
    layers :                int
                            number of hops downstream from the starting neuron
    client :                caveclient.frameworkclient.CAVEclientFull
                            CAVEclient to query information from. If not passed to the function it will create its own.

    Returns
    -------
    None

    """
    if not client:
        client = cv.CAVEclient('fanc_production_mar2021')
    def make_csvs_in_list(neuron_list):
        downstream_neurons = []
        for neuron_id in neuron_list:
            print(neuron_id, end=" ")
            if str(neuron_id)+"_downstreampartners.csv" in files:
                print("file already exists")
                pass
            else:
                downstream_new = make_csv(neuron_id, percentage_threshold, connection_threshold, folder, client)
                downstream_neurons = downstream_neurons+downstream_new
                if len(downstream_new):
                    print("created file")
                else:
                    print("no downstream partners")
        return downstream_neurons

    folder = str(start_neuron)+"-"+str(percentage_threshold)+"/"
    if not os.path.exists(folder):
        os.makedirs(folder)
    files = os.listdir(folder)

    layer_names = ["first layer", "second layer", "third layer", "fourth layer", "fifth layer"]
    if layers > 3:
        print("Warning - querying more than 3 layers can take a very long time")
    if layers > 5:
        print("limiting layers to max = 5")
        layers = 5
    
    next_layer = [start_neuron]
    for n in range(layers):
        print(layer_names[n])
        next_layer = make_csvs_in_list(next_layer)
    
def make_csv(neuron_id, percentage_threshold=0.5, connection_threshold=3, folder="", client=None):
    """Save a dataframe in .csv format of the neurons downstream from the input neuron, ordered by percentage input.
    the dataframe is sorted by high to low percent with columns -
                            bodyId: the ID of the downstream neuron
                            weight: the number of connections to the downstream neuron
                            inputs: the total number of inputs to the downstream neuron
                            percent: the percentage of input to the downstream
                                neuron represented by the connections in the input
    Parameters
    ----------
    neuron_id :             int
                            FANC neuron ID to find the downstream partners of
    percentage_threshold :  float
                            minimum percentage input required to be included in the saved csv file
    connection_threshold :  int
                            minimum number of synapses required to be considered a connection between neurons
    folder :                str
                            name of the folder in which to save the .csv file.
    client :                caveclient.frameworkclient.CAVEclientFull
                            CAVEclient to query information from. If not passed to the function it will create its own.

    Returns
    -------
    list
                        List of the IDs of downstream neurons contained in the dataframe.

    """
    table = downstream_of(neuron_id,connection_threshold)
    dataframe = table.loc[table.percent >= percentage_threshold]
    if dataframe.size == 0:
        return []
    dataframe = dataframe.rename(columns={'count': 'weight'})
    dataframe.index.names = ['bodyId']
    dataframe.to_csv(folder+str(neuron_id)+"_downstreampartners.csv")
    return dataframe.index.to_list()

def update_fanc_id(fanc_id, client=None):
    """Query CAVEclient for the newest ID associated with a neuron
    Parameters
    ----------
    neuron_id :             int or str
                            FANC neuron ID to find the current ID of
    client :                caveclient.frameworkclient.CAVEclientFull
                            CAVEclient to query information from. If not passed to the function it will create its own.

    Returns
    -------
    int or None
                        Current FANC neuron ID, or None if input ID in not valid

    """
    if not client:
        client = cv.CAVEclient('fanc_production_mar2021')
    if fanc_id == "" or fanc_id == None or fanc_id == "NotAssigned":
        return None
    newest_id = client.chunkedgraph.suggest_latest_roots(fanc_id)
    return newest_id

def mf_match(manc_id, seatable=None):
    """Query the matching table on Seatable for the FANC ID associated with a MANC neuron
    Parameters
    ----------
    manc_id :       int or str
                    MANC neuron ID to find the FANC match of
    seatable :      seaserpent.base.Table
                    Seaserpent table to query information from. If not passed to the function it will create its own.

    Returns
    -------
    int or str
                Current FANC neuron ID, or String with both connectivity match and nBlast match if they aren't the same.

    """
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

def fm_match(fanc_id, seatable=None, client=None):
    """Query the matching table on Seatable for the MANC ID associated with a FANC neuron
    Parameters
    ----------
    fanc_id :       int or str
                    FANC neuron ID to find the MANC match of
    seatable :      seaserpent.base.Table
                    Seaserpent table to query information from. If not passed to the function it will create its own.
    client :        caveclient.frameworkclient.CAVEclientFull
                    CAVEclient to query information from. If not passed to the function it will create its own.

    Returns
    -------
    int or list
                Current MANC neuron ID, or List with possible matches.

    """
    if not client:
        client = cv.CAVEclient('fanc_production_mar2021')
    if not seatable:
        seatable = ss.Table(table='fanc851_manc_nblast95_60')
    # search in the manual matches
    past_ids = [str(x) for x in client.chunkedgraph.get_past_ids(fanc_id)["past_id_map"][fanc_id]]
    past_ids.append(str(fanc_id))
    rows = seatable[seatable.manualAssignment.isin(past_ids)]
    if rows.shape[0] == 1:
        print("manually assigned")
        return int(rows.iloc[0].queryID)

    # if there are no manual matches look for nBlast and connectivity matches instead
    if rows.shape[0] == 0:
        rows = seatable[seatable.nBlastMatchID.isin(past_ids) or seatable.conMatchID.isin(past_ids)]
    
    # if there's one row and the nblast and connection matches are the same
    if rows.shape[0] == 1 and rows.iloc[0].nBlastMatchID == rows.iloc[0].conMatchID:
        return int(rows.iloc[0].queryID)
    else:
        print("no definite match for", fanc_id, ", options returned as a list.")
    return rows.queryID.to_list()

def rank_matches(fanc_csv, manc_csv, seatable=None):
    """Create a table from one MANC neuron table and one FANC table for comparison
    Parameters
    ----------
    fanc_csv :      pandas DataFrame or str
                    either a DataFrame with FANC IDs in the leftmost column or a filepath
                        (in string format) to an equivalent .csv file
    manc_csv :      pandas DataFrame or str
                    either a DataFrame with MANC IDs in the leftmost column or a filepath
                        (in string format) to an equivalent .csv file
    seatable :      seaserpent.base.Table
                    Seaserpent table to query information from. If not passed to the function it will create its own.

    Returns
    -------
    pandas DataFrame
                Combined table with all MANC IDs and the FANC IDs that could be matched to them

    """
    if not seatable:
        seatable = ss.Table(table='fanc851_manc_nblast95_60')
    if isinstance(fanc_csv, str):
        fanc = pd.read_csv(fanc_csv, index_col=0)
    else:
        fanc = fanc_csv
    if isinstance(manc_csv, str):
        manc = pd.read_csv(manc_csv, index_col=0)
    else:
        manc = manc_csv

    assigned_fanc = []
    clean_fanc = []
    for id, _ in manc.iterrows():
        fanc_id = mf_match(id,seatable)
        assigned_fanc.append(fanc_id)
        if not isinstance(fanc_id, tuple):
            clean_fanc.append(fanc_id)
        else:
            clean_fanc.append(None)
    manc["fanc"] = clean_fanc
    fanc.index = fanc.index.astype(float)
    combined_table = manc.merge(fanc, left_on="fanc", right_index=True, how="outer")
    return combined_table
