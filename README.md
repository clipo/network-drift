# network-drift

Cultural transmission models in python using simuPOP, SciPy, and NumPy. 

network-drift is a python-based software layer on top of simuPOP by Bo Peng, currently version 1.9. See http://simupop.sourceforge.net/ for source code, license, and documentation.

#Directory Structure

demography : python modules for the population structure. Currently, this consists just of network configuration. Networks can be configured as GML files or by specifiying numbers of nodes and connectivity (producing small-world graphs).

models :  python script to run the simulation models.

simuPOP_examples : python scripts that demonstrate features of simuPOP from documentation and cookbooks. 

testdata : sample GML files

#Runtime Dependency

#Module Dependencies

find . -name "*.py" | xargs grep -h 'import ' | grep -v simu | sort | cut -d' ' -f2 | uniq > required-modules.txt

#Author

Carl P. Lipo and Mark E. Madsen Copyright 2018-2019. All rights reserved. This software is made available under the Apache Software License (see file LICENSE), which allows you to use the software for commercial or non-commercial purposes, but you must attribute authorship, and derivatives must allow the user to find the original code and license.
