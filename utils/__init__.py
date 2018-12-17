# Copyright (c) $today.year.  Mark E. Madsen <mark@madsenlab.org>
#
# This work is licensed under the terms of the Creative Commons-GNU General Public Llicense 2.0, as "non-commercial/sharealike".  You may use, modify, and distribute this software for non-commercial purposes, and you must distribute any modifications under the same license.  
#
# For detailed license terms, see:
# http://creativecommons.org/licenses/GPL/2.0/

from utils.utils import init_count_traits_in_subpops
from utils.utils import constructUniformAllelicDistribution
from utils.utils import classify
from utils.utils import check_k_and_migration_rates
from utils.utils import count_traits_in_subpops
from utils.utils import init_acumulators
from utils.utils import update_acumulator
from utils.utils import update_richness_acumulator
from utils.utils import calculateAlleleAndGenotypeFrequencies
from utils.utils import mean_confidence_interval
from utils.utils import setup_output
from utils.utils import save_parameters

import math

__author__ = 'mark'

# Function for testing the partial or total ordering of a list of numbers

def strictly_increasing(L):
    return all(x<y for x, y in zip(L, L[1:]))

def strictly_decreasing(L):
    return all(x>y for x, y in zip(L, L[1:]))

def non_increasing(L):
    return all(x>=y for x, y in zip(L, L[1:]))

def non_decreasing(L):
    return all(x<=y for x, y in zip(L, L[1:]))


def simulation_burnin_time(popsize, innovrate):
    """
    Calculates burnin time, and rounds it to the nearest 1000 generation interval.

    :param popsize:
    :param innovrate:
    :return:
    """
    tmp = (9.2 * popsize) / (innovrate + 1.0) # this is conservative given the original constant is for the diploid process
    return int(math.ceil(tmp / 1000.0)) * 1000