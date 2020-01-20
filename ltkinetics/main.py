"""ltkinetics is a package for simulating nitrogenase kinetics

The differential equations were adapted from the doctoral thesis of Phillip E.
Wilson, Brigham Young University, 2005. 

written by:
Zachary Mathe, doctoral student
DeBeer Group, Max Planck Institute for Chemical Energy Conversion
10 January 2020
"""


from time import time
from copy import deepcopy
import numpy as np
from scipy.integrate import odeint


class NitrogenaseRxn:
    """for simulating a nitrogenase reaction"""
    def __init__(self):
        """Define classical LT ks, then overwrite any provided as ks_custom"""
        self.ks = {
            'k1': 5e7, 'km1': 15, 
            'k2': 200, 
            'k3': 4.4e6, 'km3': 6.4, 
            'k4': 3e6, 
            'k6': 1.2e9, 'km6': 1.75, 
            'k7': 200, # 200 in the 1984 papers, 250 in the Wilson thesis
            'k8': 8, 
            'k9': 400, 
            'k10': 4e5 , 
            'km10': 8e4, 
            'k11': 2.2e6, 'km11': 3e6, 
            }

        # Define kH constants with which to calculate gas concentrations
        self.set_temp(23)

        # Define all 62 species names
        self.y_names = [ 
            'N2', 'H2', 'NH3', 'S2O4', 'SO2', 'HSO3', 
            'E0', 'E0Fi', 'E0DF0', 'E0TF1', 'E0eDF0', 
            'HE1', 'HE1Fi', 'HE1DF0', 'HE1TF1', 'HE1eDF0', 
            'HE2', 'HE2Fi', 'HE2DF0', 'HE2TF1', 'HE2eDF0', 
            'HE3', 'HE3Fi', 'HE3DF0', 'HE3TF1', 'HE3eDF0', 
            'HE4', 'HE4Fi', 'HE4DF0', 
            'NE3', 'NE3Fi', 'NE3DF0', 'NE3TF1', 'NE3eDF0', 
            'NE34', 'NE34Fi', 'NE34DF0', 'NE34TF1', 'NE34eDF0', 
            'NE4', 'NE4Fi', 'NE4DF0', 'NE4TF1', 'NE4eDF0', 
            'NE5', 'NE5Fi', 'NE5DF0', 'NE5TF1', 'NE5eDF0', 
            'NE6', 'NE6Fi', 'NE6DF0', 'NE6TF1', 'NE6eDF0', 
            'NE7', 'NE7Fi', 'NE7DF0', 'NE7TF1', 'NE7eDF0', 
            'DF0', 'TF1', 'Fi'
            ]


    def set_ks(self, ks_custom):
        """Change kinetic constants to custom values."""
        self.ks = {**self.ks, **ks_custom}


    def set_temp(self, temp):
        """Calculate Henry's Law constants self.kH_N2 and self.kH_H2."""
        Tk = temp + 273.15
        Tstar = Tk / 100
        rho_h2o = 0.99821 - 0.000256*(Tk - 293.15)

        A_N2 = -67.3877
        B_N2 = 86.3213
        C_N2 = 24.7981
        X_N2 = np.exp(A_N2 + B_N2/Tstar + C_N2*np.log(Tstar))
        self.kH_N2 = X_N2 / (1 - X_N2) * rho_h2o * 1000 / 18.0152

        A_H2 = -48.1611
        B_H2 = 55.2845
        C_H2 = 16.8893
        X_H2 = np.exp(A_H2 + B_H2/Tstar + C_H2*np.log(Tstar))
        self.kH_H2 = X_H2/(1 - X_H2) * rho_h2o * 1000/18.0152


    def setup_y0(self, initials):
        """Setup y0 given units of atm and M and FePactive ∈ (0,1)."""
        self.initials = initials
        self.y0 = {name: 0.0 for name in self.y_names}
        default_initials = {
            'MoFe': 100e-6,
            'FeP': 100e-6,
            'DT': 10e-3,
            'FePactive': 0.45,
            'P_N2': 0.0,
            'P_H2': 0.0,
            }
        initials = {**default_initials, **initials}

        # Define gas concentrations using Henry's Law.
        self.y0['N2'] = self.kH_N2 * initials['P_N2']
        self.y0['H2'] = self.kH_H2 * initials['P_H2']

        # Define initial dithionite dissociation equilibrium.
        k6 = self.ks['k6']
        km6 = self.ks['km6']
        DT = initials['DT']
        self.y0['S2O4'] = (DT + (km6/k6 - np.sqrt((km6/k6)**2
                                 + 4*km6/k6*DT))/2)
        self.y0['SO2'] = (-1*km6/k6 + np.sqrt((km6/k6)**2 + 4*km6/k6*DT))

        # Define initial iron protein concentrations.
        self.y0['E0'] = initials['MoFe']
        self.y0['TF1'] = initials['FeP'] * initials['FePactive']
        self.y0['Fi'] = initials['FeP'] * (1 - initials['FePactive'])


    def setup_t(self, t):
        """Set self.t to provided time array."""
        self.t = t


    def integrate(self,
                  odeint_args={'atol':1e-10, 'rtol':1e-10},
                  verbose=True
                  ):
        """Integrate with odeint, having run setup_y0() and setyp_t() first.

        As of January 2020, the SciPy docs recommend using solve_ivp 
        instead of odeint. However, I have found that odeint is up to 200x 
        faster, despite the fact that both functions seem to call the same
        ODEPACK Fortran solver. 
        """
        time0 = time()
        self.sol = odeint(self.ltmodel,
                          list(self.y0.values()),
                          self.t,
                          **odeint_args
                          ).transpose()
        if verbose:
            clock = time() - time0
            print('**integrated in '+'{:.3f}'.format(clock)+' seconds**')
        self.unpack_sol(self.sol)


    def ltmodel(self, y, t):
        """
        Return dy, the derivative of all species concentrations
        t = time as a float; y = all species as a list of floats
        """
        # Unpack y and k and initialize derivatives to zero.
        [
            N2, H2, NH3, S2O4, SO2, HSO3,
            E0, E0Fi, E0DF0, E0TF1, E0eDF0,
            HE1, HE1Fi, HE1DF0, HE1TF1, HE1eDF0,
            HE2, HE2Fi, HE2DF0, HE2TF1, HE2eDF0,
            HE3, HE3Fi, HE3DF0, HE3TF1, HE3eDF0,
            HE4, HE4Fi, HE4DF0,
            NE3, NE3Fi, NE3DF0, NE3TF1, NE3eDF0,
            NE34, NE34Fi, NE34DF0, NE34TF1, NE34eDF0,
            NE4, NE4Fi, NE4DF0, NE4TF1, NE4eDF0,
            NE5, NE5Fi, NE5DF0, NE5TF1, NE5eDF0,
            NE6, NE6Fi, NE6DF0, NE6TF1, NE6eDF0,
            NE7, NE7Fi, NE7DF0, NE7TF1, NE7eDF0,
            DF0, TF1, Fi,
            ] = y

        k1 = self.ks['k1']
        km1 = self.ks['km1']
        k2 = self.ks['k2']
        k3 = self.ks['k3']
        km3 = self.ks['km3']
        k4 = self.ks['k4']
        k6 = self.ks['k6']
        km6 = self.ks['km6']
        k7 = self.ks['k7']
        k8 = self.ks['k8']
        k9 = self.ks['k9']
        k10 = self.ks['k10']
        km10 = self.ks['km10']
        k11 = self.ks['k11']
        km11 = self.ks['km11']

        [
            d_N2, d_H2, d_NH3, d_S2O4, d_SO2, d_HSO3,
            d_E0, d_E0Fi, d_E0DF0, d_E0TF1, d_E0eDF0,
            d_HE1, d_HE1Fi, d_HE1DF0, d_HE1TF1, d_HE1eDF0,
            d_HE2, d_HE2Fi, d_HE2DF0, d_HE2TF1, d_HE2eDF0,
            d_HE3, d_HE3Fi, d_HE3DF0, d_HE3TF1, d_HE3eDF0,
            d_HE4, d_HE4Fi, d_HE4DF0,
            d_NE3, d_NE3Fi, d_NE3DF0, d_NE3TF1, d_NE3eDF0,
            d_NE34, d_NE34Fi, d_NE34DF0, d_NE34TF1, d_NE34eDF0,
            d_NE4, d_NE4Fi, d_NE4DF0, d_NE4TF1, d_NE4eDF0,
            d_NE5, d_NE5Fi, d_NE5DF0, d_NE5TF1, d_NE5eDF0,
            d_NE6, d_NE6Fi, d_NE6DF0, d_NE6TF1, d_NE6eDF0,
            d_NE7, d_NE7Fi, d_NE7DF0, d_NE7TF1, d_NE7eDF0,
            d_DF0, d_TF1, d_Fi
            ] = [0.0]*62

        # Define derivatives
        d_S2O4 = k6*SO2**2 - km6*S2O4
        d_HSO3 = k4*DF0*SO2
        d_SO2 = -2*d_S2O4 - d_HSO3

        d_H2 = k7*HE2 + k8*HE3 + k9*HE4 + k10*N2*HE3 + k11*N2*HE4
        d_NH3 = km3*(NE6eDF0 + NE7eDF0)

        d_E0Fi = k3*E0*Fi - km3*E0Fi
        d_E0DF0 = k3*E0*DF0 - km3*E0DF0
        d_E0TF1 = k1*E0*TF1 - (km1 + k2)*E0TF1
        d_E0eDF0 = k2*E0TF1 - km3*E0eDF0

        d_HE1 = (km3*E0eDF0 + km3*(HE1Fi + HE1DF0) - k3*HE1*(Fi + DF0) 
                 + km1*HE1TF1 - k1*HE1*TF1
                 + k8*HE3)
        d_HE1Fi = k3*HE1*Fi - km3*HE1Fi
        d_HE1DF0 = k3*HE1*DF0 - km3*HE1DF0
        d_HE1TF1 = k1*HE1*TF1 - (km1 + k2)*HE1TF1
        d_HE1eDF0 = k2*HE1TF1 - km3*HE1eDF0

        d_HE2 = (km3*HE1eDF0 + km3*(HE2Fi + HE2DF0) - k3*HE2*(Fi + DF0)
                 + km1*HE2TF1 - k1*HE2*TF1
                 + k9*HE4
                 - k7*HE2)
        d_HE2Fi = k3*HE2*Fi - km3*HE2Fi
        d_HE2DF0 = k3*HE2*DF0 - km3*HE2DF0
        d_HE2TF1 = k1*HE2*TF1 - (km1 + k2)*HE2TF1
        d_HE2eDF0 = k2*HE2TF1 - km3*HE2eDF0

        d_HE3 = (km3*HE2eDF0 + km3*(HE3Fi + HE3DF0) - k3*HE3*(Fi + DF0)
                 + km1*HE3TF1 - k1*HE3*TF1
                 - k8*HE3 
                 - k10*N2*HE3 + km10*H2*NE3)
        d_HE3Fi = k3*HE3*Fi - km3*HE3Fi
        d_HE3DF0 = k3*HE3*DF0 - km3*HE3DF0
        d_HE3TF1 = k1*HE3*TF1 - (km1 + k2)*HE3TF1
        d_HE3eDF0 = k2*HE3TF1 - km3*HE3eDF0

        d_HE4 = (km3*HE3eDF0 + km3*(HE4Fi + HE4DF0) - k3*HE4*(Fi + DF0)
                 - k9*HE4
                 - k11*N2*HE4 + km11*H2*NE4)
        d_HE4Fi = k3*HE4*Fi - km3*HE4Fi
        d_HE4DF0 = k3*HE4*DF0 - km3*HE4DF0

        d_NE3 = (km3*(NE3Fi + NE3DF0) - k3*NE3*(Fi + DF0)
                 + km1*NE3TF1 - k1*NE3*TF1
                 + k10*N2*HE3 - km10*H2*NE3)
        d_NE3Fi = k3*NE3*Fi - km3*NE3Fi
        d_NE3DF0 = k3*NE3*DF0 - km3*NE3DF0
        d_NE3TF1 = k1*NE3*TF1 - (km1 + k2)*NE3TF1
        d_NE3eDF0 = k2*NE3TF1 - km3*NE3eDF0

        d_NE34 = (km3*NE3eDF0 + km3*(NE34Fi + NE34DF0) - k3*NE34*(Fi + DF0)
                  + km1*NE34TF1 - k1*NE34*TF1)
        d_NE34Fi = k3*NE34*Fi - km3*NE34Fi
        d_NE34DF0 = k3*NE34*DF0 - km3*NE34DF0
        d_NE34TF1 = k1*NE34*TF1 - (km1 + k2)*NE34TF1
        d_NE34eDF0 = k2*NE34TF1 - km3*NE34eDF0

        d_NE4 = (km3*(NE4Fi + NE4DF0) - k3*NE4*(Fi + DF0)
                 + km1*NE4TF1 - k1*NE4*TF1
                 + k11*N2*HE4 - km11*H2*NE4)
        d_NE4Fi = k3*NE4*Fi - km3*NE4Fi
        d_NE4DF0 = k3*NE4*DF0 - km3*NE4DF0
        d_NE4TF1 = k1*NE4*TF1 - (km1 + k2)*NE4TF1
        d_NE4eDF0 = k2*NE4TF1 - km3*NE4eDF0

        d_NE5 = (km3*(NE34eDF0 + NE4eDF0) + km3*(NE5Fi + NE5DF0) 
                 - k3*NE5*(Fi + DF0)
                 + km1*NE5TF1 - k1*NE5*TF1)
        d_NE5Fi = k3*NE5*Fi - km3*NE5Fi
        d_NE5DF0 = k3*NE5*DF0 - km3*NE5DF0
        d_NE5TF1 = k1*NE5*TF1 - (km1 + k2)*NE5TF1
        d_NE5eDF0 = k2*NE5TF1 - km3*NE5eDF0

        d_NE6 = (km3*NE5eDF0 + km3*(NE6Fi + NE6DF0) - k3*NE6*(Fi + DF0)
                 + km1*NE6TF1 - k1*NE6*TF1)
        d_NE6Fi = k3*NE6*Fi - km3*NE6Fi
        d_NE6DF0 = k3*NE6*DF0 - km3*NE6DF0
        d_NE6TF1 = k1*NE6*TF1 - (km1 + k2)*NE6TF1
        d_NE6eDF0 = k2*NE6TF1 - km3*NE6eDF0

        d_NE7 = (km3*NE6eDF0 + km3*(NE7Fi + NE7DF0) - k3*NE7*(Fi + DF0)
                 + km1*NE7TF1 - k1*NE7*TF1)
        d_NE7Fi = k3*NE7*Fi - km3*NE7Fi
        d_NE7DF0 = k3*NE7*DF0 - km3*NE7DF0
        d_NE7TF1 = k1*NE7*TF1 - (km1 + k2)*NE7TF1
        d_NE7eDF0 = k2*NE7TF1 - km3*NE7eDF0

        d_TF1 = (k4*DF0*SO2
                 + km1*(E0TF1 + HE1TF1 + HE2TF1 + HE3TF1 + NE3TF1 + NE34TF1
                        + NE4TF1 + NE5TF1 + NE6TF1 + NE7TF1)
                 - k1*TF1*(E0 + HE1 + HE2 + HE3 + NE3 + NE34
                           + NE4 + NE5 + NE6 + NE7))

        d_E0 = -1*(d_E0Fi + d_E0DF0 + d_E0TF1 + d_E0eDF0
                   + d_HE1 + d_HE1Fi + d_HE1DF0 + d_HE1TF1 + d_HE1eDF0
                   + d_HE2 + d_HE2Fi + d_HE2DF0 + d_HE2TF1 + d_HE2eDF0
                   + d_HE3 + d_HE3Fi + d_HE3DF0 + d_HE3TF1 + d_HE3eDF0
                   + d_HE4 + d_HE4Fi + d_HE4DF0 
                   + d_NE3 + d_NE3Fi + d_NE3DF0 + d_NE3TF1 + d_NE3eDF0
                   + d_NE34 + d_NE34Fi + d_NE34DF0 + d_NE34TF1 + d_NE34eDF0
                   + d_NE4 + d_NE4Fi + d_NE4DF0 + d_NE4TF1 + d_NE4eDF0
                   + d_NE5 + d_NE5Fi + d_NE5DF0 + d_NE5TF1 + d_NE5eDF0
                   + d_NE6 + d_NE6Fi + d_NE6DF0 + d_NE6TF1 + d_NE6eDF0
                   + d_NE7 + d_NE7Fi + d_NE7DF0 + d_NE7TF1 + d_NE7eDF0)

        d_Fi = -1*(d_E0Fi + d_HE1Fi + d_HE2Fi + d_HE3Fi + d_HE4Fi
                   + d_NE3Fi + d_NE34Fi + d_NE4Fi + d_NE5Fi + d_NE6Fi + d_NE7Fi)

        d_DF0 = -1*(d_TF1 + d_E0DF0 + d_E0TF1 + d_E0eDF0
                    + d_HE1DF0 + d_HE1TF1 + d_HE1eDF0
                    + d_HE2DF0 + d_HE2TF1 + d_HE2eDF0
                    + d_HE3DF0 + d_HE3TF1 + d_HE3eDF0
                    + d_HE4DF0
                    + d_NE3DF0 + d_NE3TF1 + d_NE3eDF0
                    + d_NE34DF0 + d_NE34TF1 + d_NE34eDF0
                    + d_NE4DF0 + d_NE4TF1 + d_NE4eDF0
                    + d_NE5DF0 + d_NE5TF1 + d_NE5eDF0
                    + d_NE6DF0 + d_NE6TF1 + d_NE6eDF0
                    + d_NE7DF0 + d_NE7TF1 + d_NE7eDF0)

        # Repackage and return dy
        dy = [
            d_N2, d_H2, d_NH3, d_S2O4, d_SO2, d_HSO3,
            d_E0, d_E0Fi, d_E0DF0, d_E0TF1, d_E0eDF0,
            d_HE1, d_HE1Fi, d_HE1DF0, d_HE1TF1, d_HE1eDF0,
            d_HE2, d_HE2Fi, d_HE2DF0, d_HE2TF1, d_HE2eDF0,
            d_HE3, d_HE3Fi, d_HE3DF0, d_HE3TF1, d_HE3eDF0,
            d_HE4, d_HE4Fi, d_HE4DF0,
            d_NE3, d_NE3Fi, d_NE3DF0, d_NE3TF1, d_NE3eDF0,
            d_NE34, d_NE34Fi, d_NE34DF0, d_NE34TF1, d_NE34eDF0,
            d_NE4, d_NE4Fi, d_NE4DF0, d_NE4TF1, d_NE4eDF0,
            d_NE5, d_NE5Fi, d_NE5DF0, d_NE5TF1, d_NE5eDF0,
            d_NE6, d_NE6Fi, d_NE6DF0, d_NE6TF1, d_NE6eDF0,
            d_NE7, d_NE7Fi, d_NE7DF0, d_NE7TF1, d_NE7eDF0,
            d_DF0, d_TF1, d_Fi
            ]
        return dy


    def unpack_sol(self, sol):
        """Return all species, composite E states and acid quench products

        'Enoe' is a dictionary of binned E states for which MoFe:FeP complexes
        E(n)eDF0 have been excluded. These are species in which the electron
        transfers have all occurred, but, according to older literature, the 
        accompanying proton has not yet arrived at FeMoco). 

        'E' is a dictionary of binned E states, for which the MoFe:FeP
        complexes E(n)eDF0 are binned into the subsequent E(n+1) states.
        """
        y = dict(zip(self.y_names, sol))
        self.y = y

        Enoe = {}
        Enoe['E0'] = y['E0'] + y['E0Fi'] + y['E0DF0'] + y['E0TF1']
        Enoe['E1'] = y['HE1'] + y['HE1Fi'] + y['HE1DF0'] + y['HE1TF1']
        Enoe['E2'] = y['HE2'] + y['HE2Fi'] + y['HE2DF0'] + y['HE2TF1']
        Enoe['E3'] = y['HE3'] + y['HE3Fi'] + y['HE3DF0'] + y['HE3TF1']
        Enoe['E4'] = y['HE4'] + y['HE4Fi'] + y['HE4DF0']
        Enoe['NE3'] = y['NE3'] + y['NE3Fi'] + y['NE3DF0'] + y['NE3TF1']
        Enoe['NE34'] = y['NE34'] + y['NE34Fi'] + y['NE34DF0'] + y['NE34TF1']
        Enoe['NE4'] = y['NE4'] + y['NE4Fi'] + y['NE4DF0'] + y['NE4TF1']
        Enoe['NE5'] = y['NE5'] + y['NE5Fi'] + y['NE5DF0'] + y['NE5TF1']
        Enoe['NE6'] = y['NE6'] + y['NE6Fi'] + y['NE6DF0'] + y['NE6TF1']
        Enoe['NE7'] = y['NE7'] + y['NE7Fi'] + y['NE7DF0'] + y['NE7TF1']
        self.Enoe = Enoe

        E = deepcopy(Enoe)
        # E['E0']
        E['E1'] += y['E0eDF0']
        E['E2'] += y['HE1eDF0']
        E['E3'] += y['HE2eDF0']
        E['E4'] += y['HE3eDF0']
        # E['NE3'] += 
        E['NE34'] += y['NE3eDF0']
        # E['NE4'] += 
        E['NE5'] += y['NE34eDF0'] + y['NE34eDF0']
        E['NE6'] += y['NE5eDF0']
        E['NE7'] += y['NE6eDF0']
        self.E = E

        AQ = {}
        AQ['H2'] = (y['HE2'] + y['HE2Fi'] + y['HE2DF0']
                    + y['HE2TF1'] + y['HE2eDF0']
                    + y['HE3'] + y['HE3Fi'] + y['HE3DF0']
                    + y['HE3TF1'] + y['HE3eDF0']
                    + 2*(y['HE4'] + y['HE4Fi'] + y['HE4DF0']))
        AQ['NH3'] = (2*(y['NE5'] + y['NE5Fi'] + y['NE5DF0']
                        + y['NE5TF1'] + y['NE5eDF0']
                        + y['NE6'] + y['NE6Fi'] + y['NE6DF0']
                        + y['NE6TF1'] + y['NE6eDF0'])
                     + y['NE7'] + y['NE7Fi'] + y['NE7DF0']
                     + y['NE7TF1'] + y['NE7eDF0'])
        AQ['N2H4'] = (y['NE34'] + y['NE34Fi'] + y['NE34DF0']
                      + y['NE34TF1'] + y['NE34eDF0']
                      + y['NE4'] + y['NE4Fi'] + y['NE4DF0']
                      + y['NE4TF1'] + y['NE34eDF0'])
        self.AQ = AQ
