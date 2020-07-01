from __future__ import division
from collections import defaultdict, OrderedDict
from copy import deepcopy
import psutil
import simuOpt
simuOpt.setOptions(alleleType='long', optimized=True, quiet=False)
simuOpt.setOptions(numThreads=psutil.cpu_count())
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
from time import time
import sys
import os
from collections import defaultdict
import csv

global config, sim_id, script, cores

'''
Example use of script. 
parameter-sweep-for-localization.py 
     --experiment paramsweep5 --networkfile smallworld --numloci 1 --maxinittraits 100 
    --popsize 5000 --migrationfraction 0.0005 0.001 0.0025 0.005 
    --innovrate 0.00 0.001 0.005 0.01 
    --k_values 2 25 50 100 125 
    --rewiringprob 0.001 --sub_pops 200 --maxalleles 10000 --simlength 2005 --reps 5
    
example using gml file
parameter-sweep-for-localization.py
--experiment rn-sweep 
--networkfile /Users/clipo/Documents/PycharmProjects/network-drift/data/rapa_nui/ahu.gml
 --numloci 1 --innovrate 0.00 --maxinittraits 100 
 --popsize 5000 --migrationfraction 0.000001 0.000005 0.00001 0.00005 
 --k_values 5 10 25 50 75 100 125 --rewiringprob 0.0 
 --sub_pops 150 --maxalleles 10000 --simlength 2005 --reps 5
'''

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
    parser.add_argument("--reps", help="Replicated populations per parameter set", type=int, default=3)
    parser.add_argument("--networkfile", help="Name of GML file representing the  network model for this simulation",
                        required=True, type=str, default="smallworld")
    parser.add_argument("--numloci", help="Number of loci per individual (use with care)", type=int, required=True, default=1)
    parser.add_argument("--maxinittraits", help="Max initial number of traits per locus for initialization", type=int,
                        required=True, default=50)
    parser.add_argument("--innovrate", nargs='+', help="Rate(s) at which innovations occur in population as a per-locus rate", type=float, default=[])
    parser.add_argument("--simlength", help="Time at which simulation and sampling end, defaults to 3000 generations",
                        type=int, default="20")
    parser.add_argument("--popsize", help="Initial size of population for each community in the model", type=int, required=True)
    parser.add_argument("--migrationfraction", nargs='+', help="Fraction of population that migrates each time step", type=float, required=True, default=[])
    parser.add_argument("--seed", type=int, help="Seed for random generators to ensure replicability")
    parser.add_argument( "--k_values", nargs='+', type=int, help="list of k-values to explore [e.g., 2 4 20 24]", default=[])
    parser.add_argument("--sub_pops", nargs="+", help="Number of sub populations", required=True, default=[])
    parser.add_argument("--maxalleles", type=int, help="Maximum number of alleles", default=50)
    parser.add_argument("--save_figs", type=bool, help="Save figures or not?", default=False)
    parser.add_argument("--burnintime", type=int, help="How long to wait before making measurements? ", default=2000)
    parser.add_argument("--rewiringprob", type=float, help="Probability of random rewiring", default=0)

    config = parser.parse_args()

    # setup output directories for writing
    output_path = utils.setup_output(config.experiment)

    # check the k and migration rate combinations
    check = utils.check_k_and_migration_rates(config)
    if check is not True:
        print("\nProblem(s):\t %s\n" % check)
        print("Please adjust input values for k and/or migration rate and restart.\n ")
        sys.exit()
    else:
        print("\nChecked on the migration and k values -- all looks good!\n")

    # save parameters
    utils.save_parameters(str(sys.argv), config, output_path)

    # set up the frequencies for the alleles in each loci. Here assuming a uniform distribution as a starting point
    distribution = utils.constructUniformAllelicDistribution(config.maxinittraits)

    # prepare file for output
    output_data_file_name = "%s/%s-rare-trait-output.csv" % (output_path, config.experiment)
    with open(output_data_file_name, mode='w') as output_file:
        output_writer = csv.writer(output_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        output_writer.writerow(["Iteration", "k", "NumSubPops", "Migration", "InnovationRate", "Ones_Mean",
                                "Ones_Lower", "Ones_Upper", "Twos_Mean", "Twos_Lower", "Twos_Upper", "Richness_Mean",
                                "Richness_Lower", "Richness_Upper","Fst_Mean","Fst_Lower","Fst_Upper"])
        output_file.flush()
        subpop_run_values = config.sub_pops
        k_run_values = config.k_values
        mig_run_values = config.migrationfraction
        innov_run_values = config.innovrate

        iteration=-1

        for subpop in subpop_run_values:
            if k_run_values == [0]:
                k_run_values = [2, int(float(subpop) * .1), int(float(subpop) * .2),
                            int(float(subpop) * .5), int(float(subpop) * .8),
                            int(float(subpop) * .9),
                            int(subpop) - 1]

            for k in k_run_values:
                for mig in mig_run_values:
                    for innov in innov_run_values:
                        ## let us know whats happening
                        iteration += 1
                        print("Now running with subpops: %s k-value: %s mig rate: %4f innov rate: %4f" % (subpop,k,mig,innov))
                        ## these are lists of things that simuPop will do at different stages
                        init_ops = OrderedDict()
                        pre_ops = OrderedDict()
                        post_ops = OrderedDict()

                        # Construct a demographic model
                        #networkmodel = NetworkModel( networkmodel="/Users/clipo/Documents/PycharmProjects/RapaNuiSim/notebooks/test_graph.gml",
                        networkmodel = network.NetworkModel( networkmodel=config.networkfile,
                                                         simulation_id=config.experiment,
                                                         sim_length=config.simlength,
                                                         burn_in_time=config.burnintime,
                                                         initial_subpop_size=config.popsize,
                                                         migrationfraction=mig,
                                                         sub_pops=subpop,
                                                         connectedness=k, # if 0, then distance decay
                                                         save_figs=config.save_figs,
                                                         network_iteration=iteration)

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
                        init_ops['subpop_counts'] = sp.PyOperator(utils.init_count_traits_in_subpops)
                        init_ops['Sex'] = sp.InitSex()
                        init_ops['Freq'] = sp.InitGenotype(loci=list(range(config.numloci)),freq=distribution)

                        post_ops['Innovate'] = sp.KAlleleMutator(k=config.maxalleles, rates=innov, loci=sp.ALL_AVAIL)
                        post_ops['mig']=sp.Migrator(rate=networkmodel.get_migration_matrix()) #, reps=[3])
                        post_ops['Stat-fst'] = sp.Stat(structure=sp.ALL_AVAIL)
                        post_ops['Stat-richness']=sp.Stat(alleleFreq=[0], haploFreq=[0], vars=['alleleFreq','haploFreq','alleleNum', 'genoNum'])
                        post_ops['fst_acumulation'] = sp.PyOperator(utils.update_acumulator, param=['fst','F_st'])
                        post_ops['richness_acumulation'] = sp.PyOperator(utils.update_richness_acumulator, param=('alleleFreq', 'Freq of Alleles'))
                        post_ops['class_richness']=sp.PyOperator(utils.calculateAlleleAndGenotypeFrequencies, param=(config.popsize,config.numloci))
                        post_ops['count_traits_in_subpops'] = sp.PyOperator(utils.count_traits_in_subpops, param=(config.numloci,num_pops), subPops=sp.ALL_AVAIL)

                        mating_scheme = sp.RandomSelection()

                        ## go simuPop go! evolve your way to the future!
                        sim = sp.Simulator(pops, rep=config.reps)
                        sim.evolve(initOps=list(init_ops.values()), preOps=list(pre_ops.values()), postOps=list(post_ops.values()),
                               matingScheme=mating_scheme, gen=config.simlength)

                        count=0
                        for pop in sim.populations():
                            output[count] = deepcopy(pop.dvars())
                            count+=1

                        ones_point_in_time = []
                        twos_point_in_time = []
                        richness_point_in_time = []
                        fst_point_in_time = []

                        for n in range(config.reps):
                            list_of_ones = list(output[n].ones)
                            list_of_twos = list(output[n].twos)
                            list_of_richness = list(output[n].richness)
                            list_of_fst = list(output[n].fst)
                            ones_point_in_time.append(list_of_ones[2000])
                            twos_point_in_time.append(list_of_twos[2000])
                            richness_point_in_time.append(list_of_richness[2000])
                            fst_point_in_time.append(list_of_fst[2000])

                        (ones_ave, ones_min, ones_max) = utils.mean_confidence_interval(ones_point_in_time,
                                                                                    confidence=0.95)
                        (twos_ave, twos_min, twos_max) = utils.mean_confidence_interval(twos_point_in_time,
                                                                                    confidence=0.95)
                        (richness_ave, richness_min, richness_max) = utils.mean_confidence_interval(richness_point_in_time,
                                                                                                confidence=0.95)
                        (fst_ave, fst_min, fst_max) = utils.mean_confidence_interval(fst_point_in_time,
                                                                                 confidence=0.95)

                        output_writer.writerow([iteration,k,subpop,mig,innov,ones_ave,ones_min,ones_max,
                                            twos_ave,twos_min,twos_max,richness_ave,richness_min,richness_max,fst_ave,fst_min,fst_max])
                        output_file.flush()

if __name__ == "__main__":
    main()

