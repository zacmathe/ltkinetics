"""Example 3: plotting high-flux E1 through E4, with the
still-bound MoFe•FeP species separate

Zachary Mathe, doctoral student
DeBeer Group, Max Planck Institute for Chemical Energy Conversion
15 January 2020
"""


import numpy as np
import matplotlib.pyplot as plt

import ltkinetics as lt


# Initialize reactions from the NitrogenaseRxn class
highflux = lt.NitrogenaseRxn()

# Setup initials and time points
initials = {
    'MoFe': 100e-6 * 1.9, # Remember to multiply by # Mo per MoFe dimer
    'FeP': 500e-6,
    'DT': 50e-3,
    }
highflux.setup_y0(initials)
tmax = 4
highflux.setup_t(np.linspace(0, tmax, 100000))

# Run the integration
highflux.integrate()

"""Plot the four E states, with the species EneDF0 separate. 
Following Wilson's nomenclature, we assume that MoFe reduction requires
binding of ATP-FeP, electron transfers + hydrolysis, then protein
dissociation:
    En → EnTF1 → EneDF0 → E(n+1) + DF0
"""
rxn = highflux
t, Enoe, y = rxn.t, rxn.Enoe, rxn.y
fig = plt.figure('high flux, still-bound MoFe•FeP separate')
fig.clear()
m1 = fig.add_subplot(221)
m2 = fig.add_subplot(222)
m3 = fig.add_subplot(223)
m4 = fig.add_subplot(224)
specs = [
    np.zeros_like(rxn.t),
    Enoe['E0'], y['E0eDF0'],
    Enoe['E1'], y['HE1eDF0'],
    Enoe['E2'], y['HE2eDF0'],
    Enoe['E3'], y['HE3eDF0'],
    Enoe['E4'],
    ]
stack = [sum([specs[j] for j in range(0,i+1)]) for i in range(0,len(specs))]
labs = [0, 'E0', 'E0eDF0', 'E1', 'E1eDF0', 'E2', 'E2eDF0',
        'E3', 'E3eDF0',
        'E4']
lss = ['-', '--']
zcols = ['#5C5C5B', '#e50000', '#128BD8', '#B2C700', '#9a0eea']
for ax in [m1, m2]:
    for i, spec in enumerate(stack[:-1]):
        ax.fill_between(t, stack[i+1]*1e6, y2=stack[i]*1e6,
                        label=labs[i+1],
                        color=zcols[(i+1)//2], alpha=1-((i%2)/2),
                        linewidth=0,
                        )
for ax in [m3, m4]:
    for i, spec in enumerate(specs[1:]):
        ax.plot(t, spec*1e6, label=labs[i+1],
                color=zcols[(i+1)//2], ls=lss[(i)%2]
                )
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
