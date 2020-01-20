"""Example 2 using ltkinetics.py: test the effect of reducing k7 in low flux

Zachary Mathe, doctoral student
DeBeer Group, Max Planck Institute for Chemical Energy Conversion
15 January 2020
"""


import numpy as np
import matplotlib.pyplot as plt

import ltkinetics as lt


"""Initialize reactions from the NitrogenaseRxn class"""
classical = lt.NitrogenaseRxn()
slowerk7 = lt.NitrogenaseRxn()
rxns = [classical, slowerk7]

"""Setup initials and time points"""
initials = {
    'MoFe': 100e-6 * 1.9, # Remember to multiply by number of Mo per MoFe dimer
    'FeP': 300e-6,
    'DT': 50e-3,
    }
tmax = 3
for rxn in rxns:
    rxn.setup_y0(initials)
    rxn.setup_t(np.linspace(0, tmax, 100000))

"""Use set_ks to provide a dictionary of custom kinetic constants"""
custom_ks = {'k7': 100} # The classical/default value is 200
slowerk7.set_ks(custom_ks)

"""Run the integrations"""
for rxn in rxns:
    rxn.integrate()

"""Plot four En states, with species E(n-1)eDF0 included"""
fig = plt.figure('changing k7')
fig.clear()
m1 = fig.add_subplot(121)
m2 = fig.add_subplot(122)
labs = ['E0', 'E1', 'E2', 'E3', 'E4']
suffixes = [' (classical)', ' (k7=100)']
lss = ['-', '--']
zcols = ['#5C5C5B', '#e50000', '#128BD8', '#B2C700', '#9a0eea']
for i, rxn in enumerate(rxns):
    t, E, y = rxn.t, rxn.E, rxn.y
    specs = [E['E0'], E['E1'], E['E2'], E['E3'], E['E4']]
    for ax in [m1, m2]:
        for j, spec in enumerate(specs):
            ax.plot(rxn.t, spec*1e6, c=zcols[j], ls=lss[i],
                    label=labs[j]+suffixes[i])
for ax in [m1]:
    ax.set_title('FeMoco ' + '{:.0f}'.format(rxn.initials['MoFe']*1e6) + u' μM'
                 + ', FeP ' + '{:.0f}'.format(rxn.initials['FeP']*1e6) + u' μM'
                 )
    ax.set_xlim(0, tmax/10)
    ax.set_ylabel(u'μM')
for ax in [m2]:
    ax.legend(frameon=True, ncol=2, loc=1, fontsize=12)
    ax.set_xlim(0, tmax)
for ax in [m1, m2]:
    ax.set_xlabel('s')
fig.tight_layout()
