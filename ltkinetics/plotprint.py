"""
Functions for displaying NitrogenaseRxn results

Zachary Mathe, doctoral student
DeBeer Group, Max Planck Institute for Chemical Energy Conversion
15 January 2020
"""


import numpy as np
import matplotlib.pyplot as plt


def print_Econcs(rxn, times, species):
    """Pretty print concentrations of a list of species at times"""
    for time in times:
        print('\n--at', time, ' sec--')
        tidx = abs(rxn.t - time).argmin()
        for spec in species:
            print('{:2.0f}'.format(rxn.E[spec][tidx]*1e6), spec)


# def plot_concs(rxn, times, species, figname=''):
#     """Plot species at two different time ranges in two different styles"""
#     fig = plt.figure(figname)
#     fig.clear()
#     m1 = fig.add_subplot(221)
#     m2 = fig.add_subplot(222)
#     m3 = fig.add_subplot(223)
#     m4 = fig.add_subplot(224)

#     specs = ([np.zeros_like(rxn.t)]
#              + [rxn.
#              E.E0t, y.E0eDF0,
#              E.E1t, y.HE1eDF0,
#              E.E2t, y.HE2eDF0,
#              E.E3t, y.HE3eDF0,
#              E.E4t])
