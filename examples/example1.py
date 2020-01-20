"""Example 1 using ltkinetics.py: plotting high-flux E1 through E4


Zachary Mathe, doctoral student
DeBeer Group, Max Planck Institute for Chemical Energy Conversion
15 January 2020
"""


import numpy as np
import matplotlib.pyplot as plt

import ltkinetics as lt


"""initialize reactions from the NitrogenaseRxn class"""
highflux = lt.NitrogenaseRxn()

"""Make a dictionary of initial protein concentrations (M), gas pressures
(atm) and the amount of active iron protein relative to the total.
The following default values will be used for any omitted species:
    initials = {
        'MoFe': 100e-6,
        'FeP': 100e-6,
        'DT': 10e-3,
        'FePactive': 0.45,
        'P_N2': 0.0,
        'P_H2': 0.0,
        }
"""
initials = {
    'MoFe': 100e-6 * 1.9, # Remember to multiply by # Mo per MoFe dimer
    'FeP': 500e-6,
    'DT': 50e-3,
    }
highflux.setup_y0(initials)

"""Provide a monotonic array of time points. Numpy's linspace() is easy.""" 
tmax = 4
highflux.setup_t(np.linspace(0, tmax, 100000))

"""Run the integration, automatically followed by unpack_sol. unpack_sol
creates two attributes to the reaction: Enoe and E:
    * 'Enoe' is a dictionary of binned E states for which MoFe•FeP complexes
    E(n)eDF0 have been excluded. These are species in which the electron
    transfers have all occurred, but, according to older literature, the 
    accompanying proton has not yet arrived at FeMoco). 
    * 'E' is a dictionary of binned E states, for which the MoFe•FeP
    complexes E(n)eDF0 are binned into the subsequent E(n+1) states. 
"""
highflux.integrate()

"""Print a list of E state concentrations at a list of times."""
lt.print_Econcs(highflux, [0.05, 1], ['E0', 'E1', 'E2'])

"""Plot the four E states."""
rxn = highflux
t, E, y = rxn.t, rxn.E, rxn.y
specs = [np.zeros_like(rxn.t), E['E0'], E['E1'], E['E2'], E['E3'], E['E4']]
stack = [sum([specs[j] for j in range(0,i+1)]) for i in range(0,len(specs))]

fig = plt.figure('high flux')
fig.clear()
m1 = fig.add_subplot(221)
m2 = fig.add_subplot(222)
m3 = fig.add_subplot(223)
m4 = fig.add_subplot(224)
labs = [0, 'E0','E1', 'E2', 'E3', 'E4']
zcols = ['#5C5C5B', '#e50000', '#128BD8', '#B2C700', '#9a0eea']
for ax in [m1, m2]:
    for i, spec in enumerate(stack[:-1]):
        ax.fill_between(t, stack[i+1]*1e6, y2=stack[i]*1e6,
                        label=labs[i+1], color=zcols[i], linewidth=0)
for ax in [m3, m4]:
    for i, spec in enumerate(specs[1:]):
        ax.plot(t, spec*1e6, label=labs[i+1], color=zcols[i])
for ax in [m1, m3]:
    ax.set_xlim(0, tmax/10)
    ax.set_ylabel(u'μM')
for ax in [m2, m4]:
    ax.legend(frameon=True, ncol=1, loc=1, fontsize=12)
    ax.set_xlim(0, tmax)
for ax in [m3, m4]:
    ax.set_xlabel('s')
for ax in [m1]:
    ax.set_title('FeMoco ' + '{:.0f}'.format(rxn.initials['MoFe']*1e6) + u' μM'
                 + ', FeP ' + '{:.0f}'.format(rxn.initials['FeP']*1e6) + u' μM'
                 )
fig.tight_layout()
