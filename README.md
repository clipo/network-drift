<h2>Network-drift</h2>

Cultural transmission models in python using simuPOP, SciPy, and NumPy. 

network-drift is a python-based software layer on top of simuPOP by Bo Peng, currently version 1.9. See http://simupop.sourceforge.net/ for source code, license, and documentation.

<h2>Directory Structure</h2>

demography : python modules for the population structure. Currently, this consists just of network configuration. Networks can be configured as GML files or by specifiying numbers of nodes and connectivity (producing small-world graphs).

models :  python script to run the simulation models.

simuPOP_examples : python scripts that demonstrate features of simuPOP from documentation and cookbooks. 

testdata : sample GML files

Rcode : R code and data used to produce heatmap figures for visualizing  output of simulation. 

<h2>Module Dependencies</h2>

find . -name "*.py" | xargs grep -h 'import ' | grep -v simu | sort | cut -d' ' -f2 | uniq > required-modules.txt

<h2>Runtime parameters for Figures in Paper</h2>

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
* ./network-drift/models/parameter-sweep-for-localization.py --experiment paramsweep-10-reps --networkfile smallworld --numloci 1 --maxinittraits 100 --popsize 5000  --migrationfraction 0.0001 0.00025 0.0005 0.00075 0.001 0.0025 0.005 0.0075 0.01 0.025 0.05 --innovrate 0.00 --k_values 2 10 20 40 60 80 100 120 140 --rewiringprob 0.001 --sub_pops 200 --maxalleles 10000 --simlength 2005 --reps 10

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
