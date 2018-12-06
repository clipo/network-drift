# this example is from Antao's 2015 - Bioinformation with Python Cookbook
# Chapter 5 - Simulating population structure using island and stepping-stone models
# page 138

from collections import defaultdict, OrderedDict
from copy import deepcopy
import simuPOP as sp
from simuPOP import demography
num_loci = 10
pop_size = 50
num_gens = 101
num_pops = 10
migs = [0, 0.005, 0.01, 0.02, 0.05, 0.1]
init_ops = OrderedDict()
pre_ops = OrderedDict()
post_ops = OrderedDict()
pops = sp.Population([pop_size] * num_pops, loci=[1] *
    num_loci, infoFields=['migrate_to'])


def init_accumulators(pop, param):
    accumulators = param
    for accumulator in accumulators:
        if accumulator.endswith('_sp'):
            pop.vars()[accumulator] = defaultdict(list)
        else:
            pop.vars()[accumulator] = []
    return True


def update_accumulator(pop, param):
    accumulator, var = param
    if var.endswith('_sp'):
        for sp in range(pop.numSubPop()):
            pop.vars()[accumulator][sp].append(
                deepcopy(pop.vars(sp)[var[:-3]]))
    else:
        pop.vars()[accumulator].append(deepcopy(pop.vars()[var]))
    return True

init_ops['accumulators'] = sp.PyOperator(init_accumulators,
       param=['fst'])
init_ops['Sex'] = sp.InitSex()
init_ops['Freq'] = sp.InitGenotype(freq=[0.5, 0.5])
for i, mig in enumerate(migs):
   post_ops['mig-%d' % i] = sp.Migrator(demography.migrIslandRates(mig, num_pops),
   reps=[i])
post_ops['Stat-fst'] = sp.Stat(structure=sp.ALL_AVAIL)
post_ops['fst_accumulation'] = sp.PyOperator(update_accumulator, param=('fst', 'F_st'))

mating_scheme =  sp.RandomMating()

sim = sp.Simulator(pops, rep=len(migs))
sim.evolve(initOps=list(init_ops.values()),
           preOps=list(pre_ops.values()),
           postOps=list(post_ops.values()),
           matingScheme=mating_scheme,
           gen=num_gens)

import seaborn as sns
sns.set_style('white')
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
for i, pop in enumerate(sim.populations()):
   ax.plot(pop.dvars().fst, label='mig rate %.4f' % migs[i])
ax.legend(loc=2)
ax.set_ylabel('FST')
ax.set_xlabel('Generation')
