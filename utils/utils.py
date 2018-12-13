from __future__ import division
from collections import defaultdict, OrderedDict
from copy import deepcopy
import simuPOP as sp
from simuPOP import demography
import logging as log
import numpy as np
import scipy.stats
import math
from simuPOP import demography
import demography.network as network
import os
from collections import Counter

### some functions to store stats at each timestep.

def init_count_traits_in_subpops(pop):
    '''
    Count traits in the overall population - zero out the dictionary for each loci/allele
    combination in the loci/allele trait space
    :param pop: the population object - this is passed by simuPop in the PyOperator call.
    :return: True
    '''
    pop.vars()['pop_count']=defaultdict(int)
    return True

def classify(subpops,val):
    '''
    determine which class each value is in
    :param subpops: number of subpops
    :param val: the value to check
    :return: a list with the category
    '''
    res = []
    if val == 1:
        res.append('= 1')
    if val == 2:
        res.append('= 2')
    if val < (int(int(subpops)*.05)):
        res.append('< 5%')
    if val < (int(int(subpops)*.10)):
        res.append('< 10%')
    if val < (int(int(subpops)*.2)):
        res.append('< 20%')
    if val < (int(int(subpops)*.5)):
        res.append('< 50%')
    return res

def check_k_and_migration_rates(config):
    '''
    Check the combination of k and migration rates to see if any combination results in >1.0 total
    :param config: the configuration object with the parameters from the command line
    :return: output message (if a problem) or True (if no problem
    '''
    k_values = config.k_values
    output_message = ""
    count = 0
    migs = config.migrationfraction
    for k in k_values:
        for mig in migs:
            if float(k)*float(mig) >= 1.0:
                output_message += "k=%s * mig=%4f is greater than 1.0\n" % (k,mig)
                count+=1
    if count>0:
        return output_message
    else:
        return True


def count_traits_in_subpops(pop, param):
    '''
    Count the number of subpops in which each trait occurs (1-numSubPops)
    combination in the loci/allele trait space
    :param pop: the population object - this is passed by simuPop in the PyOperator call.
    :param param: in this case pass the # of loci
    :return: True
    '''
    (num_loci, numSubPops) = param
    sp.stat(pop, haploFreq=range(0,num_loci), vars=['haploFreq_sp', 'haploNum_sp'], subPops=sp.ALL_AVAIL)

    traits_in_subpops = defaultdict(int)

    # now count for all the subpops
    for subPop in range(0,numSubPops):
        key= list(pop.vars(subPop)['haploNum'].keys())
        #traits_n_counts = pop.vars(subPop)['haploNum'][key[0]]
        haplotype_count_map= list(pop.vars(subPop)['haploNum'][key[0]].keys())
        for loci_allele_tuple in haplotype_count_map:
            traits_in_subpops[str(loci_allele_tuple)] +=1

    pop.vars()['pop_count'] = traits_in_subpops

    vals=pop.vars()['pop_count'].values()
    ones=twos=fivepercent=tenpercent=twentypercent=fiftypercent=0
    for val in vals:
        if val == 1:
            ones += 1
        if val == 2:
            twos += 1
        if val < (int(int(numSubPops) * .05)):
            fivepercent +=1
        if val < (int(int(numSubPops) * .10)):
            tenpercent +=1
        if val < (int(int(numSubPops) * .2)):
            twentypercent +=1
        if val < (int(int(numSubPops) * .5)):
            fiftypercent +=1
    pop.vars()['ones'].append(ones)
    pop.vars()['twos'].append(twos)
    pop.vars()['fivepercent'].append(fivepercent)
    pop.vars()['tenpercent'].append(tenpercent)
    pop.vars()['twentypercent'].append(twentypercent)
    pop.vars()['fiftypercent'].append(fiftypercent)
    return True

def init_acumulators(pop, param):
    acumulators = param
    for acumulator in acumulators:
        if acumulator.endswith('_sp'):
            pop.vars()[acumulator] = defaultdict(list)
        else:
            pop.vars()[acumulator] = []
            pop.vars()['allele_frequencies'] = []
            pop.vars()['haplotype_frequencies'] = []
            pop.vars()['allele_count']=[]
            pop.vars()['richness'] = []
            pop.vars()['class_freq']=[]
            pop.vars()['class_count']=[]
            pop.vars()['ones']=[]
            pop.vars()['twos']=[]
            pop.vars()['fivepercent']=[]
            pop.vars()['tenpercent']=[]
            pop.vars()['twentypercent']=[]
            pop.vars()['fiftypercent'] = []
            #pop.vars()['fst_mean']
    return True

def update_acumulator(pop, param):
    acumulator, var = param
    #for acumulator, var in sorted(param.items()):
    log.debug("acumulator: %s var: %s" % (acumulator,var))
    if  var.endswith('_sp'):
        for sp in range(pop.numSubPop()):
            pop.vars()[acumulator][sp].append(deepcopy(pop.vars(sp)[var[:-3]]))
    else:
        pop.vars()[acumulator].append(deepcopy(pop.vars()[var]))

    return True

def update_richness_acumulator(pop, param):
    (acumulator, var) = param
    if  var.endswith('_sp'):
        for sp in range(pop.numSubPop()):
            pop.vars()[acumulator][sp].append(len(pop.dvars(sp).alleleFreq[0].values()))
    else:
        pop.vars()['haplotype_frequencies'].append(len(pop.dvars().haploFreq.values()))
        pop.vars()['allele_frequencies'].append(len(pop.dvars().alleleFreq.values()))
        pop.vars()['allele_count'].append(len(pop.dvars().alleleNum))
        #pop.vars()[acumulator].append(deepcopy(pop.vars()[var]))
    return True

def calculateAlleleAndGenotypeFrequencies(pop, param):
    (popsize, num_loci) = param

    sp.stat(pop, haploFreq = range(0, num_loci), vars=['haploFreq', 'haploNum'])
    #sim.stat(pop, alleleFreq = sim.ALL_AVAIL)

    keys = list(pop.dvars().haploFreq.keys())

    haplotype_map = pop.dvars().haploFreq[keys[0]]
    haplotype_count_map = pop.dvars().haploNum[keys[0]]
    num_classes = len(haplotype_map)

    #class_freq = {'-'.join(i[0]) : str(i[1]) for i in haplotype_map.items()}
    class_freq = dict()
    for k,v in haplotype_map.items():
        key = '-'.join(str(x) for x in k)
        class_freq[key] = v
    #log.debug("class_freq packed: %s", class_freq)

    class_count = dict()
    for k,v in haplotype_count_map.items():
        key = '-'.join(str(x) for x in k)
        class_count[key] = v
    pop.vars()['richness'].append(num_classes)
    pop.vars()['class_freq'].append(class_freq)
    pop.vars()['class_count'].append(class_count)

    return True

def constructUniformAllelicDistribution(numalleles):
    """Constructs a uniform distribution of N alleles in the form of a frequency list.
        Args:
            numalleles (int):  Number of alleles present in the initial population.
        Returns:
            (list):  Array of floats, giving the initial frequency of N alleles.

    """
    divisor = 100.0 / numalleles
    frac = divisor / 100.0
    distribution = [frac] * numalleles
    return distribution

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

def setup_output(experiment="test"):
    # create output directories
    path = os.getcwd()
    path = path + "/output/"
    try:
        os.mkdir(path)
    except OSError:
        print("Creation of the directory %s failed - it might already exist? " % path)
    else:
        print("Successfully created the directory %s " % path)
    path = path + experiment

    try:
        os.mkdir(path)
    except OSError:
        print("Creation of the directory %s failed - it might already exist? " % path)
    else:
        print("Successfully created the directory %s " % path)
    return path

def save_parameters(args, config, output_path):
    '''
    Save the arguments used to run the experiment into the output director
    :param args:   the sys.argv string that contains the raw input
    :param config: a list of all the possible configuration options
    :param output_path: a path to the output directory
    :return: True if completed
    '''
    f = open(output_path + "/parameters.txt", "w+")
    f.write("Arguments used: %s\n" % args)
    f.write("--experiment: %s\n" % config.experiment)
    f.write("--debug: %s\n" % config.debug)
    f.write("--reps: %s\n" % config.reps)
    f.write("--networkfile: %s\n"% config.networkfile)
    f.write("--debug: %s\n" % config.debug)
    f.write("--numloci: %s\n" % config.numloci)
    f.write("--maxinittraits: %s\n" % config.maxinittraits)
    f.write("--innovrate: %s\n" % config.innovrate)
    f.write("--simlength: %s\n" % config.simlength)
    f.write("--popsize: %s\n" % config.popsize)
    f.write("--migrationfraction: %s\n" % config.migrationfraction)
    f.write("--seed: %s\n" % config.seed)
    f.write("--k_values: %s\n" % config.k_values)
    f.write("--sub_pops: %s\n" % config.sub_pops)
    f.write("--maxalleles: %s\n" % config.maxalleles)
    f.write("--save_figs: %s\n" % config.save_figs)
    f.write("--burnintime: %s\n" % config.burnintime)
    f.write("--rewiringprob: %s\n" % config.rewiringprob)
    f.close()

    return True
