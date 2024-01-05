import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
import pandas as pd
import os

def prepare_graph_data(csv_name):
    df = pd.read_csv(csv_name, index_col=0)
    startneuron = int(csv_name.split("_")[0])
    neurons = [startneuron]+df.index.to_list()
    weighted_edges = []
    types = {}
    connections = []
    for index, row in df.iterrows():
        weighted_edges.append((startneuron, index, row["count"]/150))
        types[index] = row["type"]
    return neurons, weighted_edges, types



def show_graph(neuronIDs, edges, types):
    G = nx.DiGraph()
    labels = {}
    for ID in neuronIDs:
        node = G.add_node(ID)
    G.add_weighted_edges_from(edges)

    widths = nx.get_edge_attributes(G, 'weight')
    nodelist = G.nodes()
    plt.title('draw_networkx')
    pos=graphviz_layout(G, prog='dot')

    color_map = []

    for node in G:
        if node in types:
            print(types[node])
        # if the node represents an interneuron
        if node not in types:
            color_map.append("green")
        elif str(types[node]).startswith("IN"):
            color_map.append("blue")
        else:
            color_map.append("red")


    nx.draw_networkx_nodes(G,pos,
                       nodelist=nodelist,
                       node_size=10,
                       node_color = color_map)
    nx.draw_networkx_edges(G,pos,
                       edgelist = widths.keys(),
                       width=list(widths.values()),
                       arrows = False)
    #nx.draw_networkx_labels(G, pos=pos,)
                        #labels=types)
    plt.show()
    

csvs = [x for x in os.listdir() if x.endswith("_downstreampartners.csv")]
neurons = []
edges = []
types = {}
for csv in csvs:
    csv_neurons, csv_edges, csv_types = prepare_graph_data(csv)
    neurons=neurons+csv_neurons
    edges=edges+csv_edges
    types=types|csv_types

show_graph(neurons, edges, types)





# for i in range(5):
#     G.add_node("Child_%i" % i)
#     G.add_node("Grandchild_%i" % i)
#     G.add_node("Greatgrandchild_%i" % i)

#     G.add_edge("ROOT", "Child_%i" % i)
#     G.add_edge("Child_%i" % i, "Grandchild_%i" % i)
#     G.add_edge("Grandchild_%i" % i, "Greatgrandchild_%i" % i)


# same layout using matplotlib with no labels
# plt.title('draw_networkx')
# pos=graphviz_layout(G, prog='dot')
# nx.draw(G, pos, with_labels=True, arrows=False)
# plt.show()
#plt.savefig('nx_test.png')