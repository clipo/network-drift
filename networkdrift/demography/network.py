#!/usr/bin/env python
# Copyright (c) 2015.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Apache Software License, Version 2.0.  See the file LICENSE for details.

"""
Description here

"""
import networkx as nx
import numpy as np
import re
import math
import simuPOP as sim
import logging as log
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from itertools import product
import sys
from sklearn.preprocessing import normalize
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
from itertools import permutations
import operator
import math
import statistics

class NetworkModel(object):
    """
    NetworkModel implements a full "demographic model" in simuPOP terms,
    that is, it manages the size and existence of subpopulations, and the
    migration matrix between them.  The basic data for this model is derived
    by the creation of a random small-world NetworkX network in the (watts_strogatz_graph) or
    the creation of a network from a GML file.
    The network edges represent a set of subpopulations
    with unique ID's, and edges between them which are weighted.  The
    weights may be determined by any underlying model (e.g., distance,
    social interaction hierarchy, etc), but need to be interpreted here
    purely as the probability of individuals moving between two subpopulations,
    since that is the implementation of interaction.
    """

    def __init__(self,
                 networkmodel="smallworld",
                 simulation_id=None,
                 sim_length=1000,
                 burn_in_time=0,
                 initial_subpop_size=500,
                 migrationfraction=0.01,
                 sub_pops=10,
                 connectedness=2,
                 rewiring_prob=0.0,
                 save_figs=True,
                 network_iteration=1,
                 xy=None,
                 output_path="test"):
        """
        :param networkmodel: Name of GML file
        :param sim_length: Number of generations to run the simulation
        :return:
        """
        # BaseMetapopulationModel.__init__(self, numGens = sim_length, initSize = initial_size, infoFields = info_fields,
        # ops = ops, sub_pops = num_subpops, connectedness = connectedness, xy=xy rewiring_prob=rewiring_prob,
        # save_figs= boolean, iteration=number)
        self.networkmodel = networkmodel  # default is small world - else GML file location
        self.sim_length = sim_length
        self.info_fields = 'migrate_to'
        self.sim_id = simulation_id
        self.burn_in_time = burn_in_time
        self.init_subpop_size = initial_subpop_size
        self.migration_fraction = migrationfraction
        self.connectedness = connectedness  # default of 3
        self.sub_pops = sub_pops  # default of 5
        self.rewiring_prob = rewiring_prob  # default of 0.0
        self._cached_migration_matrix = None
        self.subpopulation_names = []
        self.save_figs=save_figs
        self.network_iteration=network_iteration
        self.xy=[]
        self.output_path = output_path

        # Parse the GML files and create a list of NetworkX objects
        self._parse_network_model()

        # Determine the initial population configuration
        self._calculate_initial_population_configuration()

        # prime the migration matrix
        if self.connectedness==self.sub_pops:
            self._cached_migration_matrix=self._spatialMigrRates()
            log.debug(self._cached_migration_matrix)
            self.connectedness=self.sub_pops-1
        elif ".gml" in self.networkmodel:
            self._cached_migration_matrix=self._calculate_migration_matrix_from_gml()
            log.debug(self._cached_migration_matrix)
        else:
            ## used the fixed migration function for now - which determines each edge
            ## note that k * migration rate must be < 1.0
            self._cached_migration_matrix = self._calculate_fixed_migration_matrix()
            #self._cached_migration_matrix = self._calculate_migration_matrix()

    ############### Private Initialization Methods ###############

    def _parse_network_model(self):
        """
        Given a file,  read the GML files (format: <name>.gml)
        and construct a NetworkX networkmodel from the GML file
        """
        if self.networkmodel == "smallworld":
            log.debug("Creating small world Watts-Strogatz network with %s nodes and %s connections " % (
            self.sub_pops, self.connectedness))
            k=self.connectedness
            if k == self.sub_pops:
                k=k-1
            network = nx.watts_strogatz_graph(int(self.sub_pops), k, self.rewiring_prob)
            self.pos = nx.spring_layout(network, iterations=25)
            log.debug("network nodes: %s", '|'.join(sorted(str(list(network.nodes)))))
            self.network = network
            self.xy=self._set_xy_coordinates()
        elif ".gml" in self.networkmodel:
            log.debug("Opening  GML file %s:", self.networkmodel)
            network = nx.read_gml(self.networkmodel)
            log.debug("network nodes: %s", '|'.join(sorted(list(network.nodes()))))
            self.network = self._create_network_edges_from_k_value(network)
        else:
            print("There's been a problem - we haven't created the network. Bailing out!\n")
            sys.exit()

        if self.save_figs == True:
            self.print_graph()
            self.save_graph()

    def _calculate_initial_population_configuration(self):
        # num subpops is just the number of vertices in the first graph slice.
        # first_time = min(self.times)
        network = self.network
        self.sub_pops = network.number_of_nodes()
        log.debug("Number of initial subpopulations: %s", self.sub_pops)
        log.debug("list of nodes: %s", list(network.nodes(data=True)))
        # subpoplation names - have to switch them to plain strings from unicode or simuPOP won't use them as subpop names
        self.subpopulation_names = list(network.nodes)

        log.debug("calc_init_config:  subpopulation names: %s", self.subpopulation_names)

    ############### Private Methods for Call() Interface ###############

    def _get_node_label(self, g, id):
        return g.node[id]["label"].encode('utf-8', 'ignore')

    def _get_id_for_subpop_name(self, pop, name):
        return pop.subPopByName(name)

    def _get_node_parent(self, g, id):
        return g.node[id]["parent_node"].encode('utf-8', 'ignore')

    def _get_subpop_idname_map(self, pop):
        names = pop.subPopNames()
        name_id_map = dict()
        for name in names:
            id = pop.subPopByName(name)
            name_id_map[id] = name
        return name_id_map

    def _calculate_fixed_migration_matrix(self):
        for (node1, node2, data) in self.network.edges(data=True):
            self.network.add_edge(node1, node2, weight=self.migration_fraction)
        g_mat = nx.to_numpy_matrix(self.network)
        #print("normed_matrix: ", g_mat)
        return g_mat.tolist()

    def _calculate_migration_matrix_from_gml(self):
        g_mat = nx.to_numpy_matrix(self.network).astype(np.float)
        return g_mat.tolist()

    def _calculate_migration_matrix(self):
        g_mat = nx.to_numpy_matrix(self.network, weight=self.migration_fraction).astype(np.float)
        #print("g_mat: ", g_mat)
        # get the column totals
        rtot = np.sum(g_mat, axis=1)
        scaled = (g_mat / rtot) * self.migration_fraction
        diag = np.eye(np.shape(g_mat)[0]) * (1.0 - self.migration_fraction)
        g_mat_scaled = diag + scaled
        log.debug("scaled migration matrix: %s", g_mat_scaled.tolist())
        #print("g_mat_scaled: ", g_mat_scaled)
        return g_mat_scaled.tolist()

    def _spatialMigrRates(self):
        '''
        Return a migration matrix where migration rates between two
        subpopulations vary according to Euclidean distance between them.

        xy
            A list of (x,y) location for each subpopulation.

        r
            Migrate rate between two subpopulations is exp(-r*d_ij) where
            d_ij is the Euclidean distance between subpopulations i and j.
        '''
        r=self.migration_fraction
        xy=list(self.xy)
        #print(xy)
        nSubPop = self.sub_pops
        rate = []
        for i in range(nSubPop):
            rate.append([])
            for j in range(nSubPop):
                if i == j:
                    rate[-1].append(0)
                    continue
                d_ij = math.sqrt((xy[i][0] - xy[j][0]) ** 2 + (xy[i][1] - xy[j][1]) ** 2)
                rate[-1].append(math.exp(-1 * r * d_ij))
        self._cached_migration_matrix=rate
        return rate

    def _set_xy_coordinates(self):
        radius=100
        pts=self.sub_pops
        for x, y in product(range(0, int(radius) + 1, int(360 / pts)), repeat=2):
            if x ** 2 + y ** 2 <= radius ** 2:
                yield from set(((x, y), (x, -y), (-x, y), (-x, -y),))

    def _create_network_edges_from_k_value(self,network):
        new_network=nx.Graph()
        nearest_neighbor_distance = []
        number = 0
        node_pos = {}
        node_name = {}
        ## first build a list of locations and names for the new nodes
        for (num, data) in network.nodes(data=True):
            node_pos[num] = data['pos']
            node_name[num] = data['name']

        for num, xy in list(node_pos.items()):
            new_network.add_node(num, x=xy[0], y=xy[1], pos=(xy[0], xy[1]), name=node_name[num])

        self.pos = nx.get_node_attributes(new_network, 'pos')

        for num, xy in list(node_pos.items()):
            x1, y1 = xy
            # print("working on node %s" % num)
            node_distances = []
            sorted_node_distances = []
            # now iterate through to find the n closes networks
            ccount = 0
            for (num2, xy2) in list(node_pos.items()):
                x2, y2 = xy2
                if num != num2:
                    distance = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
                    node_distances.append(distance)

            sorted_node_distances.append(min(node_distances))
            smallest_distance = min(sorted_node_distances)
            nearest_neighbor_distance.append(smallest_distance)

        mean_nearest_neighbor_distance = statistics.mean(nearest_neighbor_distance)

        # now add the connections based on the degree of k specified
        for num,(x,y) in list(node_pos.items()):
            # print("working on node %s" % num)
            node_distances = {}
            # now iterate through to find the n closes networks
            ccount = 0
            for (num2, (x2,y2)) in list(node_pos.items()):
                if num != num2:
                    node_distances[num2] = (math.sqrt((x - x2) ** 2 + (y - y2) ** 2))
            sorted_node_distances = sorted(node_distances.items(), key=operator.itemgetter(1))
            list_of_edges_to_add = list(range(0, self.connectedness))
            current_k = self.connectedness
            # note: we dont want to add edges if they are already there
            for e in list_of_edges_to_add:
                ne, dist = sorted_node_distances[e]
                # network.add_edge(node_locations[num],node_locations[ne], weight=dist )
                weight=(dist / mean_nearest_neighbor_distance) * self.migration_fraction
                new_network.add_edge(num, ne, weight=weight)

        return new_network

    ###################### Public API #####################

    def get_info_fields(self):
        return self.info_fields

    def get_connectedness(self):
        return self.connectedness

    def get_initial_size(self):
        return [self.init_subpop_size] * self.sub_pops

    def get_subpopulation_names(self):
        return str(list(self.subpopulation_names))

    def get_subpopulation_sizes(self):
        return self.subpop_sizes

    def get_subpopulation_number(self):
        return len(self.subpopulation_names)

    def get_migration_matrix(self):
        return self._cached_migration_matrix

    def print_graph(self):
        """
        Show us the graph of the network.
        :return: nothing - should be a matplotlib plot
        """
        plt.subplot(111)
        nx.draw(self.network,self.pos,node_color=list(dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS).keys())[self.network_iteration], with_labels=True, font_weight='bold')
        plt.show()

    def save_graph(self):
        """
        Save the graph of the network.
        :return: nothing - should be saved file
        """
        name = "%s/k-%s.png" % (self.output_path,self.connectedness)
        nx.draw(self.network,self.pos,node_color=list(dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS).keys())[self.network_iteration], with_labels=True, font_weight='bold')
        plt.savefig(name)

    def __call__(self, pop):
        """
        Main public interface to this demography model.  When the model object is called in every time step,
        this method creates a new migration matrix.

        After migration, the stat function is called to inventory the subpopulation sizes, which are then
        returned since they're handed to the RandomSelection mating operator.

        If a new network slice is not active, the migration matrix from the previous step is applied again,
        and the new subpopulation sizes are returns to the RandomSelection mating operator as before.

        :return: A list of the subpopulation sizes for each subpopulation
        """
        if 'gen' not in pop.vars():
            gen = 0
        else:
            gen = pop.dvars().gen

        ######### Do the per tick processing ##########

        log.debug("========= Processing network  =============")
        # self._dbg_slice_pop_start(pop,gen)

        # update the migration matrix
        self._cached_migration_matrix = self._calculate_migration_matrix(gen)

        sim.migrate(pop, self._cached_migration_matrix)
        sim.stat(pop, popSize=True)
        # cache the new subpopulation names and sizes for debug and logging purposes
        # before returning them to the calling function
        self.subpopulation_names = sorted(str(list(pop.subPopNames())))
        self.subpop_sizes = pop.subPopSizes()
        #print(self.subpop_sizes)
        return pop.subPopSizes()

