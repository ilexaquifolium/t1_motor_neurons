import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
import pandas as pd
import os
import netgraph
import neuprint

def scatter(csv_name):
    if isinstance(csv_name, str):
        df = pd.read_csv(csv_name)
    else:
        df = csv_name
    manc_percent = df.percent_x.fillna(0).to_list()
    fanc_percent = df.percent_y.fillna(0).to_list()
    try:
        labels = df.fanc.to_list()
    except:
        labels = df.index.to_list()
    plt.scatter(manc_percent, fanc_percent)
    # add axes labels
    plt.xlabel('Percentage input of Neuron DNp13 in MANC R')
    plt.ylabel('Percentage input of Neuron DNp37 in FANC R')
    # add labels to all points
    for (mi, fi, label) in zip(manc_percent, fanc_percent, labels):
        if fi > 0 and mi > 0:
            plt.text(mi, fi, int(label), va='bottom', ha='center')
    plt.show()

def get_manc_types(neuronlist):
    client = neuprint.Client('neuprint-pre.janelia.org', dataset='vnc')
    neurons, _ = neuprint.fetch_neurons(neuronlist, client=client)
    neurons = neurons[["bodyId", "type", "somaSide", "rootSide"]].copy()
    neurons = neurons.set_index("bodyId")
    types = {}
    for ID in neuronlist:
        side = neurons.loc[ID].somaSide or neurons.loc[ID].rootSide
        typ = neurons.loc[ID].type + " " + str(side)
        types[ID] = typ
    return types

def get_fanc_types(neuronlist):
    match_table = pd.read_csv("MF_match_with_additions.csv")
    manc_ids = {}
    # make dictionary of fanc_id:manc_id
    for fanc_id in neuronlist:
        potential_manc_ids = match_table.loc[match_table["match"]==fanc_id]["queryID"].to_list()
        if len(potential_manc_ids):
            manc_ids[fanc_id] = potential_manc_ids[0]
    if not manc_ids:
        manc_types = {}
    else:
        manc_types = get_manc_types(list(manc_ids.values()))
    # make dictionary of fanc_id:manc_type
    fanc_types = {}
    for fanc_id in neuronlist:
        if fanc_id in manc_ids:
            typ = manc_types[manc_ids[fanc_id]]
        else:
            typ = "No type"
        fanc_types[fanc_id] = typ
    return fanc_types

# def get_fanc_types(neuronlist):
#     #searches the T1 leg and neck motor neuron csv files for any of the IDs
#     neck_MNs = pd.read_csv("neck_motor_neuron_table_v0.csv")
#     leg_MNs = pd.read_csv("FANC-MANC-matching - MNs.csv")
#     types = {}
#     for ID in neuronlist:
#         neck_id = neck_MNs[neck_MNs.pt_root_id == ID]
#         if len(neck_id):
#             types[ID] = str(ID)+"\nNeck MN"
#         else:
#             leg_id = leg_MNs[leg_MNs.latest_fanc_id == ID]
#             if len(leg_id):
#                 types[ID] = str(ID)+"\n"+str(leg_id.type.iloc[0])
#             else:# not found in either csv
#                 types[ID] = ""
#     return types

def prepare_graph_data(csv_name):
    df = pd.read_csv(csv_name, index_col=0)
    startneuron = int(csv_name.removesuffix("_downstreampartners.csv").split("/")[-1])
    neurons = [startneuron]+df.index.to_list()
    weighted_edges = []
    if len(str(startneuron)) > 12:
        types = get_fanc_types(neurons)
    else:
        types = get_manc_types(neurons)
    for index, row in df.iterrows():
        weighted_edges.append((startneuron, index, row["weight"]/50))
    return neurons, weighted_edges, types

def contract_edges(edges, u, v):
    # replace v with u
    edges = [(u,y,w) if x==v else (x,u,w) if y==v else (x,y,w) for x,y,w in edges]
    # find edges containing u
    contracting_edges = [(x,y,w) for x,y,w in edges if x==u or y==u]
    unaffected_edges = [(x,y,w) for x,y,w in edges if x!=u and y!=u]
    weighted_edges = []
    for n, edge in enumerate(contracting_edges):
        x,y,w = edge
        following_edges = [(a,b,w) for a,b,w in contracting_edges[n+1:] if x==a and y==b]
        if following_edges:
            following_edge=following_edges[0]
            w = w+following_edge[2]
        weighted_edges.append((x,y,w))
    dictionary = {(a,b):c for a,b,c in reversed(weighted_edges)}
    contracted_edges = [(a,b,c) for ((a,b),c) in dictionary.items()]
    return unaffected_edges+contracted_edges

def show_graph(neuronIDs, edges, types):
    G = nx.DiGraph()
    labels = {}
    for ID in neuronIDs:
        node = G.add_node(ID)
        labels[ID] = str(ID)+"\n"+str(types[ID])
    ###contract node groups###
    types_to_collapse = []#["IN00A021"]
    for typ in types_to_collapse:
        ids_to_join = [x for x,y in types.items() if y.split()[0]==typ]
        for id in ids_to_join[1:]:
            u = ids_to_join[0]
            v = id
            edges = contract_edges(edges,u,v)
            G.remove_node(v)
            del types[v]
            del labels[v]
            labels[u] = types[u]+"("+str(len(ids_to_join))+")"
    ######
    G.add_weighted_edges_from(edges)

    widths = nx.get_edge_attributes(G, 'weight')
    nodelist = G.nodes()

    colors= {}
    for node in G:
        if node not in types:
            colors[node] = "green"
        elif str(types[node]) == "Neck MN":
            # if the node represents an interneuron
            colors[node] = "black"
        elif "MN" in str(types[node]):
            print(node, types[node])
            # if the node represents an Motor Neuron
            colors[node] = "red"
        else:
            colors[node] = "blue"
    
    # ###remove leaves if they arent MNs###
    # remove = [node for node,degree in dict(G.degree()).items() if degree==1]
    # remove = [node for node in remove if colors[node]=="blue"]
    # G.remove_nodes_from(remove)
    # for node in remove:
    #     del labels[node],
    #     del colors[node],
    # ######

    # plot potition of nodes
    layout_type = "dot"# hierarchical layout
    #layout_type = "fdp"# spring layout
    pos=graphviz_layout(G, prog=layout_type)
    # remap positions to values between 0,1
    max_x = max(pos.values())[0]
    max_y = max([y for _,y in pos.values()])
    interactive_pos = {key: (2*x/max_x,y/max_y) for key,(x,y) in zip(pos.keys(),pos.values())}

    plt.figure(figsize=(6,6))
    plot_instance = netgraph.InteractiveGraph(
        G,
        prettify=False,
        scale=(1,1),
        arrows=True,
        node_size=0.5,
        node_layout=interactive_pos,
        node_labels=labels,
        node_color=colors,
        node_edge_color=colors,
        #node_label_offset=(0,.005),
        node_label_fontdict={"size":6},
        edge_width=nx.get_edge_attributes(G, "weight"),
        edge_color="black",
    )
    # To access the new node positions:
    node_positions = plot_instance.node_positions
    plt.show()
    
def operator_function(folder):
    # take all downstream partner .csv files in folder and plot them as a graph
    csvs = [x for x in os.listdir(folder) if x.endswith("_downstreampartners.csv")]
    neurons = []
    edges = []
    types = {}
    for csv in csvs:
        csv_neurons, csv_edges, csv_types = prepare_graph_data(folder+"/"+csv)
        neurons=list(set(neurons+csv_neurons))
        edges=edges+csv_edges
        types=types|csv_types
    show_graph(neurons, edges, types)
