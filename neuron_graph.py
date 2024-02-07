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
    client = neuprint.Client('neuprint.janelia.org', dataset='manc:v1.0')
    neurons, _ = neuprint.fetch_neurons(neuronlist)
    neurons = neurons[["bodyId", "type"]].copy()
    neurons = neurons.set_index("bodyId")
    types = {}
    for ID in neuronlist:
        typ = neurons.loc[ID].type
        types[ID] = str(ID)+"\n"+str(typ)
    return types


def get_fanc_types(neuronlist):
    #searches the T1 leg and neck motor neuron csv files for any of the IDs
    neck_MNs = pd.read_csv("neck_motor_neuron_table_v0.csv")
    leg_MNs = pd.read_csv("FANC-MANC-matching - MNs.csv")
    types = {}
    for ID in neuronlist:
        neck_id = neck_MNs[neck_MNs.pt_root_id == ID]
        if len(neck_id):
            types[ID] = "Neck MN"
        else:
            leg_id = leg_MNs[leg_MNs.latest_fanc_id == ID]
            if len(leg_id):
                types[ID] = leg_id.type.iloc[0]
            else:# not found in either csv
                types[ID] = " "
    return types


def prepare_graph_data(csv_name):
    df = pd.read_csv(csv_name, index_col=0)
    startneuron = int(csv_name.split("_")[0])
    neurons = [startneuron]+df.index.to_list()
    weighted_edges = []
    if len(str(startneuron)) > 12:
        types = get_fanc_types(neurons)
    else:
        types = get_manc_types(neurons)
    for index, row in df.iterrows():
        weighted_edges.append((startneuron, index, row["weight"]/150))
    return neurons, weighted_edges, types

def show_graph(neuronIDs, edges, types):
    G = nx.DiGraph()
    labels = {}
    for ID in neuronIDs:
        node = G.add_node(ID)
    G.add_weighted_edges_from(edges)

    widths = nx.get_edge_attributes(G, 'weight')
    nodelist = G.nodes()
    
    # plot potition of nodes
    pos=graphviz_layout(G, prog='dot')
    # remap positions to values between 0,1
    max_x = max(pos.values())[0]
    max_y = max([y for _,y in pos.values()])
    interactive_pos = {key: (3*x/max_x,y/max_y) for key,(x,y) in zip(pos.keys(),pos.values())}

    colors= {}
    for node in G:
        if node not in types:
            colors[node] = "green"
        elif str(types[node]) == "Neck MN":
            # if the node represents an interneuron
            colors[node] = "black"
        elif "MN" in str(types[node]):
            print(node, types[node].split("\n")[1])
            # if the node represents an Motor Neuron
            colors[node] = "red"
        else:
            colors[node] = "blue"


    # nx.draw_networkx_nodes(G,pos,
    #                    nodelist=nodelist,
    #                    node_size=10,
    #                    node_color = color_map)
    # nx.draw_networkx_edges(G,pos,
    #                    alpha=0.5,
    #                    edgelist = widths.keys(),
    #                    width=list(widths.values()),
    #                    arrows = False)
    # higher_pos = {}
    # for id in pos:
    #     x,y = pos[id]
    #     higher_pos[id] = x,y+10
    # nx.draw_networkx_labels(G, pos=higher_pos, font_size=6,
    #                     labels=types)
    plot_instance = netgraph.InteractiveGraph(
        G,
        node_size=1,
        node_layout=interactive_pos,
        node_labels=types,
        node_color=colors,
        node_edge_color=colors,
        node_label_offset=(0,.005)
    )
    # To access the new node positions:
    node_positions = plot_instance.node_positions
    plt.show()
    
def operator_function():
    csvs = [x for x in os.listdir() if x.endswith("_downstreampartners.csv")]
    male_csvs = [x for x in csvs if len(x.split("_")[0]) < 12]
    female_csvs = [x for x in csvs if len(x.split("_")[0]) > 12]
    neurons = []
    edges = []
    types = {}
    for csv in male_csvs:
        csv_neurons, csv_edges, csv_types = prepare_graph_data(csv)
        neurons=list(set(neurons+csv_neurons))
        edges=edges+csv_edges
        types=types|csv_types
    show_graph(neurons, edges, types)
