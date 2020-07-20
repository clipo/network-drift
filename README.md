<h2>Network-drift</h2>

Cultural transmission models in python using simuPOP, SciPy, and NumPy. 

network-drift is a python-based software layer on top of simuPOP by Bo Peng, currently version 1.9. See http://simupop.sourceforge.net/ for source code, license, and documentation.

<h2>Overview</h2>

The network-drift project examines the effects of drift on variability within populations of varying size and structure. 
By definition, drift (whether among genetic variants or cultural ones) is a random process. In large populations, drift
 tends to cause relatively small changes in trait frequency. In small populations, however, drift can lead to rapid changes in trait frequencies, 
 potentially resulting in fixation or loss of variants. In these situations, drift can produce dramatic differences in historic outcomes: 
 two populations drawn from a single population may begin with an identical frequency of traits but can quickly diverge in composition.  
 In addition, drift in small populations often leads to rapid loss of diversity and richness in variation.
 
 The effect of drift on the frequency of traits in a population can be easily simulated. For example, compare drift over 
 time in a population of varying size (Figs. 1 and 2).  In each simulation, traits are initially present in 50% of individuals.  
 For each population size, we track the changing frequency of traits when replication is purely subject to sampling 
 (i.e., no selection or other biases). The smaller the population, the greater the chance that traits will be eliminated. 
 Of interest here are the relative rates of change for traits going to zero. While any trait in any population size has a 
 chance of being reduced to zero over time, smaller populations tend to lose variants more rapidly than larger populations. 
 This is the basis of Kimura and Ohta’s (1969) work that showed that the time until a trait goes extinct in a population 
 depends on two parameters, Ne (the effective population) and p (the initial trait frequency).
 
 ![Drift]("https://github.com/clipo/network-drift/blob/master/images/drift-with-varying-population-sizes.png?raw=true "Drift among populations of varying sizes")
 Figure 1. Drift under varying population sizes: (A) 5000 (B) 500 (C) 100 (D) 50 (E) 25. For each population size, 100 replicates 
 are shown where traits begin at 50% and due to the effects of drift stochastically change in frequency over 1000 time steps. 
 In some replicates, traits go to fixation and in some, they go extinct.

![Drift Effects]("https://github.com/clipo/network-drift/blob/master/images/drift-effects-varying-population-sizes.png?raw=true "Drift among populations of varying sizes")
Figure 2. The effects of drift on traits in populations of 5 different sizes (5000, 500, 100, 50, 25).  This figure shows the 
relative rate of trait loss across the 1000 replicate simulations over 1000 time steps for each population size. In each replicate, 
traits begin at 50% prevalence in the population. The smaller the population size, the more traits are likely to be lost from the 
population. In the large population, 
however, no traits are lost.

One significant factor that can contribute to the impact of drift on variability is population structure. The impact of drift 
is greatest when populations are well-mixed. The greater the degree of structure within a population, the more likely that 
variability will be retained in a population. We can demonstrate the relation of population structure and its impact on 
drift by modeling population interaction as a network. Structure within a population can be represented by a network where 
vertices represent individuals (N) and edges represent the potential interaction between those individuals (e.g., mating or 
social learning). The structure of the network then varies by the number of 
edges between individual vertices (k, the network degree), from immediate neighbors (k=2) to all other vertices (k=N-1).   

![Network Structure]("https://github.com/clipo/network-drift/blob/master/images/schneider-et-al-2016-figure.png?raw=true "Drift among populations of varying sizes")

Schneider and colleagues (2016) have shown that the effect of drift on diversity can be countered by a 
combination of mutation (or innovation) rate and/or highly structured (low k) networks (Figure 7). 
Following Schneider et al. (2016), given a population of size N and mutation rate μ, 
drift dominates whenever 2μN⪡1. For 2μN⪢1, on the other hand, mutation dominates over drift, maximizing genetic diversity.
 The transition occurs at a threshold, μc = 1/2N, where the equilibrium distribution of traits frequencies becomes uniform.  
 This threshold (kc) is a critical point: above it drift is insensitive to population spatial structure. Below kc,
  small degrees of connectivity (k) afford high degrees of spatial structure when combined with mutation rate leads to
   increases in diversity (de Aguiar et al., 2009; Martins et al., 2013). As shown in Figure 6, diversity can be increased 
   either by increasing the mutation rate at fixed k or by decreasing k at fixed mutation rate. Schneider et al. (2016:15) 
   remark that “in consonance with classical results, extreme restriction in [gene] flow is required for structuring to have an 
   effect. In fact, the critical mutation rate 
above which drift is overcome changes significantly only when the degree of the network becomes very small.”

In terms of information retention within a population, our interest is less in the effect of mutation (or innovation) 
but rather on population structure (k).  While innovation would increase overall diversity, its effect would be to 
alter information and thus the potential contribution to future fitness consequences. Instead, we need to focus
 on population structure since it is conceivable that the way individuals or sub-populations within a population 
 interact could serve to retain information richness as well as diversity. In this way, we would expect populations 
 for which the retention of information has strong selective advantages to be highly structured.  To explore this possibility, 
 this simulation explores various impacts of varying configurations of populations. To accommodate the aggregate nature of the 
 archaeological record, these models need to accommodate populations interacting at two different scales: individuals interacting 
 within a small 
community with nominal structure and individual communities interacting differentially amongst one another.  
		
We use simuPOP (Peng and Kimmel 2005) to simulate drift within populations of varying structure. SimuPOP is a python-based population 
simulator that allows one to evolve populations forward in time under varying configurations of mutation, recombination, 
migration, and population/subpopulation sizes. We based our simulations on a simple Wright–Fisher model (Fisher 1923; Wright 1931) 
that explores changes in a haploid population of fixed size, N. Traits are modeled as values within a single loci, in a way 
that is equivalent to attributes along a single dimension (sensu Dunnell 1971). In the model, traits for individuals are 
derived each time-step by sampling with replacement from the pool of other individuals (i.e., the previous “generation”). 
This pattern of copying traits is effectively random, implying that an individual has equal probability to interact with a
nyone else in the population. At each step in which an individual copies traits, there is a fixed chance of innovation in 
which a new trait is introduced. 

To accommodate population structure, we divided the overall population into a number of subpopulations. Within each subpopulation, 
copying is assumed to be random, but is not allowed between subpopulations. With subpopulations, drift can produce a unique 
combination of traits even if the initial conditions begin identically. We treat copying between subpopulations as migration, 
the exchange of information and traits between subpopulations. Subpopulations can be configured to interact with other subpopulations 
by copying traits depending on a network configuration. In simuPOP, we can also vary the migration rate, the percentage 
of individuals from subpopulations that copy from other subpopulations. 

<h2>References</h2>

de Aguiar, M.A.M., Baranger, M., Baptestini, E.M., Kaufman, L., Bar-Yam, Y., 2009. Global patterns of speciation and diversity. Nature 460, 384.

Dunnell, R.C., 1971. Systematics in prehistory. Free Press, New York.

Fisher, R.A., 1923. XXI.—On the Dominance Ratio. Proceedings of the Royal Society of Edinburgh 42, 321–341. https://doi.org/10.1017/S0370164600023993

Kimura, M., Ohta, T., 1969. The Average Number of Generations until Fixation of a Mutant Gene in a Finite Population. Genetics 61, 763–771.

Martins, A.B., Aguiar, M.A.M. de, Bar-Yam, Y., 2013. Evolution and stability of ring species. PNAS 110, 5080–5084. https://doi.org/10.1073/pnas.1217034110

Peng, B., Kimmel, M., 2005. simuPOP: a forward-time population genetics simulation environment. Bioinformatics 21, 3686–3687. https://doi.org/10.1093/bioinformatics/bti584

Schneider, D.M., Martins, A.B., de Aguiar, M.A.M., 2016. The mutation–drift balance in spatially structured populations. Journal of Theoretical Biology 402, 9–17. https://doi.org/10.1016/j.jtbi.2016.04.024

Wright, S., 1931. Evolution in Mendelian Populations. Genetics 16, 97–159.


<h2>Directory Structure</h2>

networkdrift/demography : python modules for the population structure. Currently, this consists just of network configuration. Networks can be configured as GML files or by specifiying numbers of nodes and connectivity (producing small-world graphs).

models :  python scripts that run the simulation models.

simuPOP_examples : python scripts that demonstrate features of simuPOP from documentation and cookbooks. 

data/rapa_nui : contains shapefiles and gml file for ahu locations on Rapa Nui

testdata : sample GML files

Rcode : R code and data used to produce heatmap figures for visualizing  output of simulation. 

images : images in this README.md document

<h2>Module Dependencies</h2>

find . -name "*.py" | xargs grep -h 'import ' | grep -v simu | sort | cut -d' ' -f2 | uniq > required-modules.txt

<h2>Scripts</h2>

There are 3 primary python scripts that use the demography modules to run simulations under varying conditions. These 
are found in the /models directory.

<h4>network-k-eval.py</h4>
This python script runs simulations of a population broken into a series of fixed subpopulations that vary in their 
interaction configuration (k). k vlaues can vary from 2 (i.e., each sub population interacts with 2 of its neighbors) to
 a sccenario in which each subpopulation interacts with all other subpopulations (N-1). This script produces 
a series of visualizations that include the network and various metrics (richness, diversity, rare trait numbers) shown over time.

Arguments taken by the python script include:
* --experiment  (Name for experiment, is used as name of output directory for figures and data that are generated.)
* --debug (turn on debugging output, default=False)
* --reps (Determines the number repeated runs that are done to produce confidence intervals. 
Number of replicated populations per parameter set", default=1)
* --networkfile (Name of GML file representing the  network model for this simulation. Use "smallworld" for generic configurations. 
In cases where one wants to simulation using geographic locations, provide a GML file that has positions for each location. [e.g., ahu.gml].))
* --numloci (Number of loci per individual. Usually 1.)
* --maxinittraits (Maximum initial number of traits per locus for initialization).
* --innovrate (Rate at which innovations occur in population as a per-locus rate, default=0.001)
* --simlength (Number of generations. Point at which simulation and sampling end)
* --popsize (Initial size of population)
* --migrationfraction (Fraction of population that interacts with another subpopulation at each time step, default=0.0001)
* --seed", type=int, help="Seed for random generators to ensure replicability")
* --k_values (list of k-values to explore [e.g., 2 4 20 24])
* --sub_pops (Number of sub populations, default=10)
* --maxalleles (Maximum number of alleles", default=50)
* --save_figs (Save figures or not?, default=True)
* --burnintime (How long to wait before making measurements? Doesn't do anything in this script context, default=200)
* --rewiringprob (Probability of random rewiring of network", default=0)

<h4>network-subpop-eval.py</h4>
This python script runs simulations of a population broken into a varying number of subpopulations that have a fixed
interaction configuration (k). Subpopulations can vary from 2 to the size of the total population. This script produces 
a series of visualizations that include the network and various metrics (richness, diversity, rare trait numbers) shown over time.  
 
Arguments taken by the python script include:
* --experiment  (Name for experiment, is used as name of output directory for figures and data that are generated.)
* --debug (turn on debugging output, default=False)
* --reps (Determines the number repeated runs that are done to produce confidence intervals. 
Number of replicated populations per parameter set", default=1)
* --networkfile (Name of GML file representing the  network model for this simulation. Use "smallworld" for generic configurations. In 
cases where one wants to simulation using geographic locations, provide a GML file that has positions for each location. [e.g., ahu.gml].)
* --numloci (Number of loci per individual. Usually 1.)
* --maxinittraits (Maximum initial number of traits per locus for initialization).
* --innovrate, (Rate at which innovations occur in population as a per-locus rate, default=0.001)
* --simlength (Number of generations. Point at which simulation and sampling end)
* --popsize (Initial size of population)
* --migrationfraction (Fraction of population that interacts with another subpopulation at each time step, default=0.0001)
* --seed", type=int, help="Seed for random generators to ensure replicability")
* --k_values (k-value to explore [e.g., 20 ])
* --sub_pops (List of sub population numbers for simulation runs. [e.g., 10 100 200])
* --maxalleles (Maximum number of alleles", default=50)
* --save_figs (Save figures or not?, default=True)
* --burnintime (How long to wait before making measurements? Doesn't do anything in this script context,, default=2000)
* --rewiringprob (Probability of random rewiring of network", default=0)

<h4>parameter-sweep-for-localization.py</h4>
This python script runs parameter sweeps of various values that include the k, number of subpops, migration rate, and innovation rate.  
The output of the script is a CSV file that shows the metrics that each configuration reached at the point of time in that is set by the "burnintime" parameter.
These CSV files then can be read by the R code (in repository) to produce heatmap figures. Note that running this script can
 take fairly long depending on the number of configuration combinations provided. 
 
Arguments taken by the python script include:
* --experiment  (Name for experiment, is used as name of output directory for figures and data that are generated.)
* --debug (turn on debugging output, default=False)
* --reps (Determines the number repeated runs that are done to produce confidence intervals. 
Number of replicated populations per parameter set", default=1)
* --networkfile (Name of GML file representing the  network model for this simulation. Use "smallworld" for generic configurations. In 
cases where one wants to simulation using geographic locations, provide a GML file that has positions for each location. [e.g., ahu.gml].)
* --numloci (Number of loci per individual. Usually 1.)
* --maxinittraits (Maximum initial number of traits per locus for initialization).
* --innovrate, (List of rates at which innovations occur in population as a per-locus rate [e.g., 0.1 0.001 0.0001])
* --simlength (Number of generations. Point at which simulation and sampling end)
* --popsize (Initial size of population)
* --migrationfraction (list of fraction of population that interacts with another subpopulation at each time step, default=0.0001)
* --seed", type=int, help="Seed for random generators to ensure replicability")
* --k_values (List of k-values to explore [e.g., 20 50 150])
* --sub_pops (List of sub population numbers for simulation runs. [e.g., 10 100 200])
* --maxalleles (Maximum number of alleles", default=50)
* --save_figs (Save figures or not?, default=True)
* --burnintime (How long to wait before making measurements? Doesn't do anything in this script context,, default=2000)
* --rewiringprob (Probability of random rewiring of network", default=0

<h2>Runtime parameters for figures included in paper</h2>

The following runtime parameters were used to generate the figures in the associated paper.

<h3>Figures 3-4</h3>
Examines the effect of varying degrees of network connectivity (k) on the richness and diversity of traits. 
* ./network-drift/models/network-k-eval.py
--experiment smallworld-k-2-50-120 --networkfile smallworld --numloci 1 --maxinittraits 100 --popsize 5000 --migrationfraction 0.0001 --innovrate 0.00 --k_values 2 50 120 --rewiringprob 0.001 --sub_pops 150 --maxalleles 10000 --simlength 2000 --reps 10 --save_figs True

<h3>Figures 5-6</h3>
Examines the effect of varying the number of sub-populations within a population on the richness and diversity of traits.
*./network-drift/models/network-subpop-eval.py
--experiment smallworld-subpops-2-25-100 --networkfile smallworld --numloci 1 --maxinittraits 100 --popsize 5000 --migrationfraction 0.0001 --innovrate 0.00 --k_values 2 --rewiringprob 0.001 --sub_pops 2 50 200 --maxalleles 10000 --simlength 2000 --reps 10 --save_figs True

<h3>Figures 7-12</h3>
Examines the effect of varying migration rates between subpopulations as well as varying degrees of connectivity (k). 
Note: this process generates the data for these figures. The R-code produces the heatmap visualizations
* ./network-drift/models/parameter-sweep-for-localization.py --experiment paramsweep-10-reps --networkfile smallworld --numloci 1 --maxinittraits 100 --popsize 5000  --migrationfraction 0.0005 0.001 0.0025 0.005 --innovrate 0.00 --k_values 2 25 50 100 125 --rewiringprob 0.001 --sub_pops 200 --maxalleles 10000 --simlength 2005 --reps 10

<h3>Figure 13 (Creates network configurations for Figure 13)</h3>
Produces the varying degrees of connectivity among ahu locations on Rapa Nui. 
* ./network-drift/models/shapefile_to_gml.py
--shapefile ../data/rapa_nui/ahu.shp --migrationfraction 0.0001 --connectedness 5 --output rapa_nui_network --nodecolor blue
* ./network-drift/models/shapefile_to_gml.py
--shapefile ../data/rapa_nui/ahu.shp --migrationfraction 0.0001 --connectedness 50 --output rapa_nui_network --nodecolor green
* ./network-drift/models/shapefile_to_gml.py
--shapefile ../data/rapa_nui/ahu.shp --migrationfraction 0.0001 --connectedness 120 --output rapa_nui_network --nodecolor red

<h3>Figures 14-15</h3>
Evaluates effects on diversity and richness of traits under varying degrees of connectivity Rapa Nui ahu locations
* ./network-drift/models/network-k-eval.py --experiment rapa_nui-k-5-50-140 --networkfile ./network-drift/data/rapa_nui/ahu.gml --numloci 1 --maxinittraits 100 --popsize 5000 --migrationfraction .0001 --innovrate 0.00 --k_values 5 50 140 --rewiringprob 0.0001 --sub_pops 150 --maxalleles 10000 --simlength 2000 --reps 10  --save_figs True

<h3>Figures 16-17</h3>
Explores varying migration rates and connectivity (k) among ahu locations on Rapa Nui. Note: this process generates the data for these figures. The R-code produces the heatmap visualizations
* ./network-drift/models/parameter-sweep-for-localization.py --experiment rn-sweep --networkfile ./network-drift/data/rapa_nui/ahu.gml --numloci 1 --innovrate 0.00 --maxinittraits 100 --popsize 5000 --migrationfraction 0.000001 0.000005 0.00001 0.00005 --k_values 5 10 25 50 75 100 125 --rewiringprob 0.0 --sub_pops 150 --maxalleles 10000 --simlength 2005 --reps 5

<h2>Authors</h2>

Carl P. Lipo and Mark E. Madsen Copyright 2018-2020. All rights reserved. This software is made available under the Apache Software License (see file LICENSE), which allows you to use the software for commercial or non-commercial purposes, but you must attribute authorship, and derivatives must allow the user to find the original code and license.
