from __future__ import division
from collections import defaultdict, OrderedDict
from copy import deepcopy
import simuOpt
simuOpt.setOptions(alleleType='long', optimized=True, quiet=False)
import math
import seaborn as sns
sns.set_style('white')
import matplotlib.pyplot as plt
import demography.network as network
import utils.utils as utils
import simuPOP as sp
from simuPOP import demography
import logging as log
import numpy as np
import scipy.stats
import argparse
import uuid
from matplotlib import colors as mcolors
from time import time

global config, sim_id, script, cores

# now set up the basic parameters of the simulation (need to change this to a config file...)
# num_loci = 4        # for now we just need one
# pop_size = 5000
# num_gens = 100
migs = [0.001, 0.01, 0.1]
# pop_list = [100, 500, 1000]
# innovation_rate = 0.005
# MAXALLELES = 10000
# connectedness=3 ## k
# sub_pops=25
# migration_fraction=0.01
# num_starting_alleles=1000
# save_state=True
# burn_in_time = 4000

output={}

k_values=[2]

def setup(parser):
    config = parser.parse_args()
    sim_id = uuid.uuid1().urn
    script = __file__

    if config.debug == '1':
        log.basicConfig(level=log.DEBUG, format='%(asctime)s %(levelname)s: %(message)s')
    else:
        log.basicConfig(level=log.INFO, format='%(asctime)s %(levelname)s: %(message)s')

def main():
    start = time()
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", help="provide name for experiment", required=True, type=str)
    parser.add_argument("--debug", help="turn on debugging output")
    parser.add_argument("--reps", help="Replicated populations per parameter set", type=int, default=1)
    parser.add_argument("--networkfile", help="Name of GML file representing the  network model for this simulation",
                        required=True, type=str)
    parser.add_argument("--numloci", help="Number of loci per individual", type=int, required=True)
    parser.add_argument("--maxinittraits", help="Max initial number of traits per locus for initialization", type=int,
                        required=True)
    parser.add_argument("--innovrate", help="Rate at which innovations occur in population as a per-locus rate", type=float, default=0.001)
    parser.add_argument("--simlength", help="Time at which simulation and sampling end, defaults to 3000 generations",
                        type=int, default="3000")
    parser.add_argument("--popsize", help="Initial size of population for each community in the model", type=int, required=True)
    parser.add_argument("--migrationfraction", help="Fraction of population that migrates each time step", type=float, required=True, default=0.2)
    parser.add_argument("--seed", type=int, help="Seed for random generators to ensure replicability")
    parser.add_argument("--k_values", type=list, help="list of k-values to explore [2,4,20,25]", required=True, default=[2])
    parser.add_argument("--sub_pops", type=int, help="Number of sub populations", required=True, default=10)
    parser.add_argument("--maxalleles", type=int, help="Maximum number of alleles", default=1000)
    parser.add_argument("--save_figs", type=bool, help="Save figures or not?", default=True)
    parser.add_argument("--burnintime", type=int, help="How long to wait before making measurements? ", default=2000)

    config = parser.parse_args()

    run_param=k_values
    ## initialize the output dictionary
    for k in run_param:
        output[k]=[]

    # set up the frequencies for the alleles in each loci. Here assuming a uniform distribution as a starting point
    distribution = utils.constructUniformAllelicDistribution(config.maxinittraits)

    iteration_number=-1
    for param_value in run_param:
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
                                             sub_pops=config.sub_pops,
                                             connectedness=param_value, # if 0, then distance decay
                                             save_figs=config.save_figs,
                                             network_iteration=iteration_number)

        num_pops = networkmodel.get_subpopulation_number()
        sub_pop_size = int(config.popsize / num_pops)

        # The regional network model defines both of these, in order to configure an initial population for evolution
        # Construct the initial population
        pops = sp.Population( size = networkmodel.get_initial_size(),
                             subPopNames = str(list(networkmodel.get_subpopulation_names())),
                             infoFields = 'migrate_to',
                             ploidy=1,
                             loci=config.numloci )

        ### now set up the activities
        init_ops['acumulators'] = sp.PyOperator(utils.init_acumulators, param=['fst','alleleFreq', 'haploFreq'])
        init_ops['Sex'] = sp.InitSex()

        init_ops['Freq'] = sp.InitGenotype(loci=list(range(config.numloci)),freq=distribution)

        post_ops['Innovate']=sp.KAlleleMutator(k=config.maxalleles, rates=config.innovrate, loci=sp.ALL_AVAIL)
        #post_ops['mig'] = sp.Migrator(demography.migrIslandRates(migration_rate, num_pops)) #, reps=[i])
        post_ops['mig']=sp.Migrator(rate=networkmodel.get_migration_matrix())
        #for i, mig in enumerate(migs):
        #        post_ops['mig-%d' % i] = sp.Migrator(demography.migrIslandRates(mig, num_pops), reps=[i])

        post_ops['Stat-fst'] =sp.Stat(structure=sp.ALL_AVAIL)
        #post_ops['haploFreq']=sp.stat(pops, haploFreq=[0], vars=['haploFreq', 'haploNum'])
        #post_ops['alleleFreq']=sp.stat(pops, alleleFreq=sp.ALL_AVAIL)

        post_ops['Stat-richness']=sp.Stat(alleleFreq=[0], haploFreq=[0], vars=['alleleFreq','haploFreq','alleleNum', 'genoNum'])
        post_ops['fst_acumulation'] = sp.PyOperator(utils.update_acumulator, param=['fst','F_st'])
        post_ops['richness_acumulation'] = sp.PyOperator(utils.update_richness_acumulator, param=('alleleFreq', 'Freq of Alleles'))
        post_ops['class_richness']=sp.PyOperator(utils.sampleAlleleAndGenotypeFrequencies, param=(config.popsize,config.numloci))

        mating_scheme = sp.RandomSelection()
        #mating_scheme=sp.RandomSelection(subPopSize=sub_pop_size)

        ## go simuPop go! evolve your way to the future!
        sim = sp.Simulator(pops) #, rep=3)

        sim.evolve(initOps=list(init_ops.values()), preOps=list(pre_ops.values()), postOps=list(post_ops.values()),
                   matingScheme=mating_scheme, gen=config.simlength)

        print("now evolving... k= %s" % param_value)

        # # now make a figure of the Fst results
        # fig = plt.figure(figsize=(16, 9))
        # ax = fig.add_subplot(111)
        # for pop, mig in zip(sim.populations(), migs):
        #     ax.plot(pop.dvars().fst, label='Migration rate %.4f' % mig)
        # ax.legend(loc=2)
        # ax.set_ylabel('FST')
        # ax.set_xlabel('Generation')
        # plt.show()
        #
        # # # now make a figure of richness.
        # fig2 = plt.figure(figsize=(16, 9))
        # ax = fig2.add_subplot(111)
        # for pop, mig in zip(sim.populations(), migs):
        #     ax.plot(pop.dvars().richness, label='Migration rate %.4f' % mig)
        # ax.legend(loc=2)
        # ax.set_ylabel('Richness')
        # ax.set_xlabel('Generation')
        # plt.show()
        # output[param_value]=deepcopy(pop.dvars())
        # #print(output)

    sum_fig = plt.figure(figsize=(16,9))
    ax=sum_fig.add_subplot(111)
    iteration=-1
    for k in run_param:
        iteration += 1
        ax.plot(output[k].fst, color=list(dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS).keys())[iteration], label='k - %s' % k)
    ax.legend(loc=2)
    ax.set_ylabel('Fst')
    ax.set_xlabel('Generations')
    plt.show()
    sum_fig.savefig('sum_fig.png', bbox_inches='tight')

    rich_fig = plt.figure(figsize=(16,9))
    ax=rich_fig.add_subplot(111)
    iteration=-1
    for k in run_param:
        iteration+=1
        ax.plot(output[k].richness, color=list(dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS).keys())[iteration], label='k - %s' % k)
    ax.legend(loc=2)
    ax.set_ylabel('Richness')
    ax.set_xlabel('Generations')
    plt.show()
    rich_fig.savefig('richness.png', bbox_inches='tight')

    ## output CI for the parameters
    if len(run_param) > 0:
        print("k value      mean            lower                upper")
        for k in run_param:
            print("k - %s: %s" % (k, utils.mean_confidence_interval(output[k].fst[4000:8000], confidence=0.95)))

# # now make a figure of the haplotypeFreq results...
# fig3 = plt.figure(figsize=(16, 9))
# ax = fig3.add_subplot(111)
# for pop, mig in zip(sim.populations(), migs):
#     ax.plot(pop.dvars().alleleNum, label='Migration rate %.4f' % mig)
# ax.legend(loc=2)
# ax.set_ylabel('Allele Numbers')
# ax.set_xlabel('Generation')
# plt.show()
# #
# # now make a figure of the haplotypeFreq results...
# fig4 = plt.figure(figsize=(16, 9))
# ax = fig4.add_subplot(111)
# for pop, mig in zip(sim.populations(), migs):
#     ax.plot(pop.dvars().class_richness, label='Migration %.4f' % mig)
# ax.legend(loc=2)
# ax.set_ylabel(' Class Richness')
# ax.set_xlabel('Generation')
# plt.show()

if __name__ == "__main__":
    main()