
import simuPOP as sim
from simuPOP.utils import migrIslandRates

p = [0.2, 0.3, 0.5]

pop = sim.Population(size=[10000]*3, loci=1, infoFields='migrate_to')
simu = sim.Simulator(pop, rep=2)

simu.evolve(
    initOps=[sim.InitSex()] +
        [sim.InitGenotype(prop=[p[i], 1-p[i]], subPops=i) for i in range(3)],

    preOps=sim.Migrator(rate=migrIslandRates(0.01, 3), reps=0),
    matingScheme=sim.RandomMating(),
    postOps=[
        sim.Stat(alleleFreq=0, structure=0, vars='alleleFreq_sp', step=50),
        sim.PyEval(r"'Fst=%.3f (%s)' % (F_st, ', '.join(['%.2f’ % "
                   "subPop[x]['alleleFreq'][0][0] for x in range(3)]))",
                   step=50),
        sim.PyOutput('\n', reps=-1, step=50),
    ],
    gen=201
)

