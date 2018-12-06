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

def update_acumulator(pop, param):
    acumulator, var = param
    #for acumulator, var in sorted(param.items()):
    log.debug("acumulator: %s var: %s" % (acumulator,var))
    if  var.endswith('_sp'):
        for sp in range(pop.numSubPop()):
            pop.vars()[acumulator][sp].append(deepcopy(pop.vars(sp)[var[:-3]]))
    else:
        pop.vars()[acumulator].append(deepcopy(pop.vars()[var]))
        #pop.vars()['richness'].append(deepcopy(pop.vars()['haploFreq']))
        #output[run_param].append(deepcopy(pop.vars()[var]))
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