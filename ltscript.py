#!/usr/bin/env python3
"""Script for running ltkinetics from the command line

Simply modify the values in the first section, then execute the script
in the terminal. 

Provide the 'MoFe' concentration multiplied by the amount of FeMoco per MoFe
dimer. In the literature, this can be 1.4 for Kp and 1.9 for Av nitrogenases.

The output file has 12 columns, starting with time:
t, E0, E1, E2, E3, E4, E3+N2, E4+N2 (from E3+N2), E4+N2 (from E4), E5, E6, E7
"""


import pprint
import numpy as np
import matplotlib.pyplot as plt

import ltkinetics as lt

pp = pprint.PrettyPrinter()


##################################
##      Requried Variables      ##
##  •concentrations in M        ##
##  •pressures in atm           ##
##  •active FeP ratio ∈ (0,1)   ##
##################################
initials = {
    'MoFe': 50e-6 * 1.9, # Nominally 50 μM, times 1.9 for vinelandii Mo/MoFe
    'FeP': 200e-6, 
    'DT': 50e-3,
    'FePactive': 0.45, 
    'P_N2': 1.0,
    'P_H2': 0.0,
    }

tmax = 5 # Simulation will run from 0 to tmax seconds

name = 'ltscript-example' # Name to prepend to output files



##################################
##      Optional Variables      ##
##  •number of time points      ##
##  •kinetic constants (SI)     ##
##################################
tpoints = 10000 

custom_ks = { 
    'k1': 5e7,
    'km1': 15, 
    'k2': 200, 
    'k3': 4.4e6,
    'km3': 6.4, 
    'k4': 3e6, 
    'k6': 1.2e9,
    'km6': 1.75, 
    'k7': 200, 
    'k8': 8, 
    'k9': 400, 
    'k10': 4e5 , 
    'km10': 8e4, 
    'k11': 2.2e6,
    'km11': 3e6, 
    }



##################################
##    Run sim & write output    ##
##################################
print('\n######## ltkinetics script run ########')
print('\n**initial conditions**')
pp.pprint(initials)
print('\n**LT kinetic constants**')
pp.pprint(custom_ks)
# print('\n')

r = lt.NitrogenaseRxn()
r.setup_y0(initials)
r.setup_t(np.linspace(0, tmax, tpoints))
r.integrate()

filename = name+'-E-pops.dat'
data = np.array(r.t)
for E in r.E.values():
    data = np.column_stack((data, E))
np.savetxt(filename, data, fmt='%4.12e')
print('\n**output has been written to**\n', filename)
print('\n**ltkinetics script run finished**')
