import gdal
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
from itertools import permutations
import operator
import math


migration_fraction=0.001
networkfile="/Users/clipo/Google Drive File Stream/Team Drives/Rapa Nui/SAA 2019/Simulation/data/ahu.shp"
network = nx.read_shp(networkfile)
k=5
## the mean nearest neighbor distance will be the scale
mean_nearest_neighbor_distance=398.6199479781432
exp_nearest_neighbor_distance=759.4643281637791

#print(list(network.nodes(data=True)))

number=0
node_locations={}
for ((x,y), data) in network.nodes(data=True):
    node_locations[number]=(x,y)
    number += 1
#    network.add_edge(node1, node2, weight=migration_fraction)

new_network=nx.Graph()
for num,xy in list(node_locations.items()):
    new_network.add_node(num, x=xy[0], y=xy[1], pos=(xy[0],xy[1]))
pos=nx.get_node_attributes(new_network,'pos')
for num,xy in list(node_locations.items()):
    #print("working on node %s" % num)
    node_distances = {}
    # now iterate through to find the n closes networks
    ccount=0
    for (num2,xy2) in list(node_locations.items()):
        if num != num2:
            node_distances[num2]=(math.sqrt((xy[0]-xy2[0])**2 + (xy[1]-xy2[1])**2))
    sorted_node_distances = sorted(node_distances.items(), key=operator.itemgetter(1))
    list_of_edges_to_add=list(range(0,k))
    current_k=k
    for e in list_of_edges_to_add:
        ne, dist = sorted_node_distances[e]
        if new_network.has_edge(num,ne):
            current_k += 1
            list_of_edges_to_add.append(current_k)
        else:
            #network.add_edge(node_locations[num],node_locations[ne], weight=dist )
            new_network.add_edge(num, ne, weight=dist/exp_nearest_neighbor_distance*migration_fraction)
        #print("adding edge from %s to %s" % (num, ne))

nx.draw(new_network, pos)
nx.write_gml(new_network,"/Users/clipo/Google Drive File Stream/Team Drives/Rapa Nui/SAA 2019/Simulation/data/rapa_nui-k5.gml")
plt.show()