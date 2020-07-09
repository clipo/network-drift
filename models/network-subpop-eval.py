from __future__ import division
from collections import defaultdict, OrderedDict
from copy import deepcopy
import simuOpt
simuOpt.setOptions(alleleType='long', optimized=True, quiet=False)
import math
import seaborn as sns
sns.set_style('white')
import matplotlib.pyplot as plt
import networkdrift.demography.network as network
from networkdrift.utils import utils
import simuPOP as sp
import logging as log
import numpy as np
import scipy.stats
import argparse
import uuid
from matplotlib import colors as mcolors
from time import time
import sys
import os
from collections import defaultdict

global config, sim_id, script, cores

output=defaultdict(dict)

def setup(parser):
    config = parser.parse_args()
    sim_id = uuid.uuid1().urn
    script = __file__

    if config.debug == '1':
        log.basicConfig(level=log.DEBUG, format='%(asctime)s %(levelname)s: %(message)s')
    else:
        log.basicConfig(level=log.INFO, format='%(asctime)s %(levelname)s: %(message)s')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", help="provide name for experiment", required=True, type=str, default="test")
    parser.add_argument("--debug", help="turn on debugging output")
    parser.add_argument("--reps", help="Replicated populations per parameter set", type=int, default=1)
    parser.add_argument("--networkfile", help="Name of GML file representing the  network model for this simulation",
                        required=True, type=str)
    parser.add_argument("--numloci", help="Number of loci per individual", type=int, required=True)
    parser.add_argument("--maxinittraits", help="Max initial number of traits per locus for initialization", type=int,
                        required=True)
    parser.add_argument("--innovrate", help="Rate at which innovations occur in population as a per-locus rate", type=float, default=0.001)
    parser.add_argument("--simlength", help="Time at which simulation and sampling end, defaults to 3000 generations",
                        type=int, default="20")
    parser.add_argument("--popsize", help="Initial size of population for each community in the model", type=int, required=True)
    parser.add_argument("--migrationfraction", nargs='+', help="Fraction of population that migrates each time step",
                        type=float, required=True, default=[])
    parser.add_argument("--seed", type=int, help="Seed for random generators to ensure replicability")
    parser.add_argument( "--k_values", nargs='+', type=int, help="list of k-values to explore [e.g., 2 4 20 24", default=[])
    parser.add_argument("--sub_pops", nargs='+', help="Number of sub populations", required=True, default=[10])
    parser.add_argument("--maxalleles", type=int, help="Maximum number of alleles", default=50)
    parser.add_argument("--save_figs", type=bool, help="Save figures or not?", default=True)
    parser.add_argument("--burnintime", type=int, help="How long to wait before making measurements? ", default=2000)
    parser.add_argument("--rewiringprob", type=float, help="Probability of random rewiring", default=0)

    config = parser.parse_args()

    # check the k and migration rate combinations
    for kvalue in config.k_values:
        if float(kvalue) * float(config.migrationfraction) >= 1.0:
            print("k=%s * mig=%4f is greater than 1.0\n" % (kvalue, config.migrationfraction))
            print("Please adjust input values for k and/or migration rate and restart.\n ")
        sys.exit()

    # setup output directories for writing
    output_path = utils.setup_output(config.experiment)

    # save parameters
    utils.save_parameters(str(sys.argv), config, output_path)

    k_run_values=config.k_values
    subpop_run_values = config.sub_pops

    ## make sure the k values are less than # of subpops and > 1
    for k in k_run_values:
        for subnum in subpop_run_values:
            if int(k) > int(subnum) or int(k) < 2:
                print("k values can not be greater than the number of sub populations. k = %s subpops = %s \n" % (k, subnum))
                sys.exit()

    ## initialize the output dictionary
    for k in k_run_values:
        for sb in subpop_run_values:
            output[k][sb]={}

    # set up the frequencies for the alleles in each loci. Here assuming a uniform distribution as a starting point
    distribution = utils.constructUniformAllelicDistribution(config.maxinittraits)

    iteration_number=-1
    for k in k_run_values:
        for subnum in subpop_run_values:
            iteration_number += 1
            ## these are lists of things that simuPop will do at different stages
            init_ops = OrderedDict()
            pre_ops = OrderedDict()
            post_ops = OrderedDict()

            # Construct a demographic model from a collection of network slices which represent a temporal network
            # of changing subpopulations and interaction strengths.  This object is Callable, and simply is handed
            # to the mating function which applies it during the copying process
            #networkmodel = NetworkModel( networkmodel="/Users/clipo/Documents/PycharmProjects/RapaNuiSim/notebooks/test_graph.gml",
            networkmodel = network.NetworkModel( networkmodel="smallworld",
                                                 simulation_id=config.experiment,
                                                 sim_length=config.simlength,
                                                 burn_in_time=config.burnintime,
                                                 initial_subpop_size=config.popsize,
                                                 migrationfraction=config.migrationfraction,
                                                 sub_pops=subnum,
                                                 connectedness=k, # if 0, then distance decay
                                                 save_figs=config.save_figs,
                                                 network_iteration=iteration_number)

            num_pops = networkmodel.get_subpopulation_number()
            sub_pop_size = int(config.popsize / num_pops)

            # The regional network model defines both of these, in order to configure an initial population for evolution
            # Construct the initial population
            pops = sp.Population(size = [sub_pop_size]*num_pops,
                                 subPopNames = str(list(networkmodel.get_subpopulation_names())),
                                 infoFields = 'migrate_to',
                                 ploidy=1,
                                 loci=config.numloci )

            ### now set up the activities
            init_ops['acumulators'] = sp.PyOperator(utils.init_acumulators, param=['fst','alleleFreq', 'haploFreq'])
            init_ops['Sex'] = sp.InitSex()

            init_ops['Freq'] = sp.InitGenotype(loci=list(range(config.numloci)),freq=distribution)

            post_ops['Innovate'] = sp.KAlleleMutator(k=config.maxalleles, rates=config.innovrate, loci=sp.ALL_AVAIL)
            post_ops['mig']=sp.Migrator(rate=networkmodel.get_migration_matrix()) #, reps=[3])
            #for i, mig in enumerate(migs):
            #        post_ops['mig-%d' % i] = sp.Migrator(demography.migrIslandRates(mig, num_pops), reps=[i])

            post_ops['Stat-fst'] = sp.Stat(structure=sp.ALL_AVAIL)

            post_ops['Stat-richness']=sp.Stat(alleleFreq=[0], haploFreq=[0], vars=['alleleFreq','haploFreq','alleleNum', 'genoNum'])
            post_ops['fst_acumulation'] = sp.PyOperator(utils.update_acumulator, param=['fst','F_st'])
            post_ops['richness_acumulation'] = sp.PyOperator(utils.update_richness_acumulator, param=('alleleFreq', 'Freq of Alleles'))
            post_ops['class_richness']=sp.PyOperator(utils.calculateAlleleAndGenotypeFrequencies, param=(config.popsize,config.numloci))

            mating_scheme = sp.RandomSelection()
            #mating_scheme=sp.RandomSelection(subPopSize=sub_pop_size)

            ## go simuPop go! evolve your way to the future!
            sim = sp.Simulator(pops, rep=config.reps)
            print("now evolving... k = %s with sub_pops = %s" % (k,subnum))
            sim.evolve(initOps=list(init_ops.values()), preOps=list(pre_ops.values()), postOps=list(post_ops.values()),
                       matingScheme=mating_scheme, gen=config.simlength)

            # now make a figure of the Fst results
            fig = plt.figure(figsize=(16, 9))
            ax = fig.add_subplot(111)
            count=0
            for pop in sim.populations():
                ax.plot(pop.dvars().fst, label='Replicate: %s' % count)
                output[k][subnum][count] = deepcopy(pop.dvars())
                count += 1
            ax.legend(loc=2, fontsize=24)
            ax.set_ylabel('FST', fontsize=24)
            ax.set_xlabel('Generations', fontsize=24)
            ax.tick_params(axis='x', labelsize=20)
            ax.tick_params(axis='y', labelsize=20)
            plt.show()

    sum_fig = plt.figure(figsize=(16,9))
    ax=sum_fig.add_subplot(111)
    iteration=-1
    for k in k_run_values:
        for subnum in subpop_run_values:
            iteration += 1
            # only label the first one
            for n in range(config.reps):
                if n==0:
                    ax.plot(output[k][subnum][n].fst, color=list(dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS).keys())[iteration], label='k = %s subpops = %s' % (k, subnum))
                else:
                    ax.plot(output[k][subnum][n].fst,
                            color=list(dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS).keys())[iteration])
    ax.legend(loc=2, fontsize=24)
    ax.set_ylabel('Fst', fontsize=24)
    ax.set_xlabel('Generations', fontsize=24)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    plt.show()
    savefilename= output_path + "/sum_fig.png"
    sum_fig.savefig(savefilename, bbox_inches='tight')

    rich_fig = plt.figure(figsize=(16,9))
    ax=rich_fig.add_subplot(111)
    iteration=-1
    for k in k_run_values:
        for sb in subpop_run_values:
            iteration+=1
            # only add a label for the first one (not all the replicates)
            for n in range(config.reps):
                if n==0:
                    ax.plot(output[k][subnum][n].richness,
                            color=list(dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS).keys())[iteration], label='k = %s subpops = %s' % (k, subnum))
                else:
                    ax.plot(output[k][subnum][n].richness,
                            color=list(dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS).keys())[iteration])
    ax.legend(loc=2, fontsize=24)
    ax.set_ylabel('Richness', fontsize=24)
    ax.set_xlabel('Generations', fontsize=24)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    plt.show()
    savefilename = output_path + "/richness.png"
    rich_fig.savefig(savefilename, bbox_inches='tight')

    ## output CI for the parameters

    summary_fig = plt.figure(figsize=(16, 9))
    ax = summary_fig.add_subplot(111)

    iteration = -1
    for k in k_run_values:
        for subnum in subpop_run_values:
            iteration += 1
            CI_average = []
            CI_min = []
            CI_max = []
            for t in range(len(output[k][subnum][0].fst)):
                point_in_time = []
                for n in range(config.reps):
                    list_of_points = list(output[k][subnum][n].fst)
                    point_in_time.append(list_of_points[t])
                (ave, min, max) = utils.mean_confidence_interval(point_in_time, confidence=0.95)
                CI_average.append(ave)
                CI_min.append(min)
                CI_max.append(max)
            ax.plot(list(CI_average), color=list(dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS).keys())[iteration],label='k = %s subpops = %s' % (k, subnum))
            ax.plot(list(CI_min), "--", color="0.5")
            ax.plot(list(CI_max), "--", color="0.5")
            ax.fill_between(list(CI_average), list(CI_max), list(CI_min), color="None", linestyle="--")
    ax.legend(loc=2, fontsize=24)
    ax.set_ylabel('Fst', fontsize=24)
    ax.set_xlabel('Generations', fontsize=24)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    plt.show()
    savefilename = output_path + "/summary-ci.png"
    summary_fig.savefig(savefilename, bbox_inches='tight')

    ## now the richness graph
    richness_sum_fig = plt.figure(figsize=(16, 9))
    ax = richness_sum_fig.add_subplot(111)

    iteration=-1
    for k in k_run_values:
        for subnum in subpop_run_values:
            iteration += 1
            CI_average = []
            CI_min = []
            CI_max = []
            for t in range(len(output[k][subnum][0].richness)):
                point_in_time = []
                for n in range(config.reps):
                    list_of_points = list(output[k][subnum][n].richness)
                    point_in_time.append(list_of_points[t])
                (ave, min, max) = utils.mean_confidence_interval(point_in_time, confidence=0.95)
                CI_average.append(ave)
                CI_min.append(min)
                CI_max.append(max)
            ax.plot(list(CI_average), color=list(dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS).keys())[iteration],label='k = %s subpops = %s' % (k, subnum))
            ax.plot(list(CI_min), "--", color="0.5")
            ax.plot(list(CI_max), "--", color="0.5")
            ax.fill_between(list(CI_average), list(CI_max), list(CI_min), color="None", linestyle="--")
    ax.legend(loc=2, fontsize=24)
    ax.set_ylabel('Richness', fontsize=24)
    ax.set_xlabel('Generation', fontsize=24)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    plt.show()
    savefilename = output_path + "/richness-ci.png"
    richness_sum_fig.savefig(savefilename, bbox_inches='tight')

if __name__ == "__main__":
    main()

