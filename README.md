# ltkinetics
Python code for simulating the kinetics of the nitrogenase enzyme system,
following the work of Thorneley and Lowe. 

Each reaction is an instance of the NitrogenaseRxn class. 
The user sets up the starting conditions and SciPy `odeint` (LSODA under the hood) numerically integrates the system forward in time. 
The oridinary differential equations (ODEs) were adapted from the Mathematica code in the doctoral thesis of Phillip E. Wilson (1) 


## Features
The original LT kinetic constants are used by default, but can be easily modified. 
The use of the package is demonstrated with two examples. Follow example 1 first. 


### Disclaimer
This is a work in progress by a non-programmer. All feedback is welcome. 


### References
1. [Wilson's Thesis, advised by Watt](https://scholarsarchive.byu.edu/etd/516/)
