"""Example 4: replicating figures from Wilson's thesis

Zachary Mathe, doctoral student
DeBeer Group, Max Planck Institute for Chemical Energy Conversion
15 January 2020
"""


import numpy as np
import matplotlib.pyplot as plt

import ltkinetics as lt


#{{{ figure B-3

# Initialize reactions from the NitrogenaseRxn class
ra = lt.NitrogenaseRxn()

# Setup initials
initials = {
    'MoFe': 34e-6 * 1.0, # Thesis didn't use a Mo/MoFe factor here
    'FeP': 133e-6,
    'DT': 10e-3,
    'P_N2': 0.0
    }
ra.setup_y0(initials)

# Setup time points and run the integration
tmax = 2 
ra.setup_t(np.linspace(0, tmax, 10000))
ra.integrate()

# Plot H2
rxn = ra
fig = plt.figure('wilson fig B-3')
fig.clear()
m1 = fig.add_subplot(111)
labs = ['H2 evolved', 'H2 acidquench', 'H2 total']
zcols = ['#5C5C5B', '#e50000', '#128BD8', '#B2C700', '#9a0eea']
t, E, y, AQ = rxn.t, rxn.E, rxn.y, rxn.AQ
specs = [y['H2'], AQ['H2'], y['H2']+AQ['H2']]
for j, spec in enumerate(specs):
    m1.plot(t, spec*1e6, c=zcols[j], label=labs[j])
m1.legend(frameon=False, fontsize=12)
m1.set_xlim(0, tmax)
m1.set_ylim(0, 140)
m1.set_xlabel('s')
m1.set_ylabel(u'μM')
fig.tight_layout()
#}}}

#{{{ figure 2-9

# Initialize reactions from the NitrogenaseRxn class
ra = lt.NitrogenaseRxn()
rb = lt.NitrogenaseRxn()
rxns = [ra, rb]

# Setup initials 
initials = {
    'MoFe': 33e-6 * 1.3, 
    'FeP': 100e-6,
    'DT': 10e-3,
    'P_N2': 1.0
    }
ra.setup_y0(initials)
initials = {
    'MoFe': 15.38e-6 * 1.3, 
    'FeP': 64e-6,
    'DT': 10e-3,
    'P_N2': 1.0
    }
rb.setup_y0(initials)

# Setup time points and run the integrations
tmax = 10
for rxn in rxns:
    rxn.setup_t(np.linspace(0, tmax, 10000))
    rxn.integrate()

# Plot NH3 (including that released upon acid quench)
fig = plt.figure('wilson fig 2-9')
fig.clear()
m1 = fig.add_subplot(111)
lss = ['-', '--']
for i, rxn in enumerate(rxns):
    t, y, AQ = rxn.t, rxn.y, rxn.AQ
    specs = [y['NH3']+AQ['NH3']]
    for j, spec in enumerate(specs):
        m1.plot(t, spec*1e6, c='k', ls=lss[i])
m1.set_xlim(0, tmax)
m1.set_ylim(0, 160)
m1.set_xlabel('s')
m1.set_ylabel(u'μM')
fig.tight_layout()
#}}}
