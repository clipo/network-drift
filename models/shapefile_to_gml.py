import gdal
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
from itertools import permutations
import operator
import math
import statistics
import argparse

global config

'''
Example use of script. 
python3 ./shapefile_to_gml.py 
    --shapefile ../data/rapa_nui/ahu.shp
    --migrationfraction 0.0001
    --connectedness 5
    --nodecolor blue
    --output ahu

'''

def print_graph(network,pos, nodecolor):
    """
    Show us the graph of the network.
    :return: nothing - should be a matplotlib plot
    """
    plt.subplot(111)
    nx.draw(network, pos,
            node_color=nodecolor,
            with_labels=False, font_weight='bold')
    plt.show()
    return True

def save_graph_plot(network,pos,outputname,connectedness,nodecolor):
    """
    Save the graph of the network.
    :return: nothing - should be saved file
    """
    name = "%s-%s.svg" % (outputname, connectedness)
    nx.draw(network, pos,
            node_color=nodecolor,
            with_labels=True, font_weight='bold')
    plt.savefig(name)
    return True

def save_gml(network,outputname,nodecolor):
    """
    Save the GML file
    :return: nothing - should be saved file
    """
    name = "%s.gml" % outputname
    nx.write_gml(network, name)
    return True

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", help="turn on debugging output")
    parser.add_argument("--shapefile", help="name of shapefile to import", required=True)
    parser.add_argument("--migrationfraction", help="base migration fraction", type=float, required=True, default=0.001)
    parser.add_argument("--connectedness", help="k value used to connect nodes in network", type=int, required=True, default=5)
    parser.add_argument("--output", help="name of output gml", required=True)
    parser.add_argument("--nodecolor", help="color of nodes", required=True, default="blue")

    config = parser.parse_args()

    ## the mean nearest neighbor distance will be used to scale the migration fraction - using the distance between nodes
    migration_fraction =config.migrationfraction
    networkfile = config.shapefile
    output = config.output
    nodecolor=config.nodecolor

    k=config.connectedness

    network = nx.read_shp(networkfile)
    number=0
    node_locations={}
    node_name={}

    ## first build a list of locations and names for the new nodes
    for ((x,y), data) in network.nodes(data=True):
        node_locations[number]=(x,y)
        node_name[number]=data['NAME']
        number += 1

    nearest_neighbor_distance=[]
    new_network=nx.Graph()
    for num,xy in list(node_locations.items()):
        new_network.add_node(num, x=xy[0], y=xy[1], pos=(xy[0],xy[1]), name=node_name[num])
    pos=nx.get_node_attributes(new_network,'pos')

    for num,xy in list(node_locations.items()):

        x1,y1=xy
        #print("working on node %s" % num)
        node_distances = []
        sorted_node_distances = []
        # now iterate through to find the n closes networks
        ccount=0
        for (num2,xy2) in list(node_locations.items()):
            x2,y2=xy2
            if num != num2:
                distance=math.sqrt((x1-x2)**2 + (y1-y2)**2)
                node_distances.append(distance)

        sorted_node_distances.append(min(node_distances))
        smallest_distance=min(sorted_node_distances)
        nearest_neighbor_distance.append(smallest_distance)

    mean_nearest_neighbor_distance=statistics.mean(nearest_neighbor_distance)

    for num, xy in list(node_locations.items()):
        # print("working on node %s" % num)
        node_distances = {}
        # now iterate through to find the n closes networks
        ccount = 0
        for (num2, xy2) in list(node_locations.items()):
            if num != num2:
                node_distances[num2] = (math.sqrt((xy[0] - xy2[0]) ** 2 + (xy[1] - xy2[1]) ** 2))
        sorted_node_distances = sorted(node_distances.items(), key=operator.itemgetter(1))
        list_of_edges_to_add=list(range(0,k))
        current_k=k
        for e in list_of_edges_to_add:
                ne, dist = sorted_node_distances[e]
                # network.add_edge(node_locations[num],node_locations[ne], weight=dist )
                weight = (dist / mean_nearest_neighbor_distance) * migration_fraction
                new_network.add_edge(num, ne, weight=weight)

    print_graph(new_network, pos,nodecolor)
    save_graph_plot(new_network,pos,output,k,nodecolor)
    save_gml(new_network,output,nodecolor)

if __name__ == "__main__":
    main()
