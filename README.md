# ltkinetics
A python package for simulating the kinetics of the nitrogenase enzyme system, following the work of Thorneley & Lowe. 

Each reaction is an instance of the NitrogenaseRxn class. 
The user sets up the starting conditions and SciPy's [odeint](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html) numerically integrates the system forward in time. 
The oridinary differential equations (ODEs) were adapted from the Mathematica code in the doctoral thesis of Phillip E. Wilson (1).


## Usage
ltkinetics can be run easily from the terminal using the provided ltscript.py, or imported as a regular python package. The latter option is required to access all features. 

#### With ltscript.py:
1. Open ltscript.py in a text editor and provide starting conditions.
2. Run the script using e.g. `./ltscript.py`. The script must remain in the same directory as the ltkinetics folder. 
3. Open the `*-E-pops.dat` file in any plotting program.

#### As a python package:
1. Run `pip install git+https://github.com/zacmathe/ltkinetics` (or put `ltkinetics` in your working directory, or add it to your `PYTHONPATH`). 
2. Open example1.py in your IDE of choice.


## Features
The use of the package is demonstrated with three examples. I recommend starting with example 1, which is more thoroughly commented.

In example 1, we integrate a high-flux reaction until the steady state is reached and plot the results both traditionally and stacked:
![ex1-result1](examples/ex1-result1.png)

In example 2, we demonstrate the funciton `set_ks`, which allows the modification of the default LT kinetic constants:
![ex2-result1](examples/ex2-result1.png)

In example 3, we do the same thing as example 1, but separate still-bound MoFe(red)•FeP(ox,ADP) species:
![ex3-result1](examples/ex3-result1.png)


## Requirements
ltkinetics was written in python 3.7 and may not be compatible with previous versions of python. Your python version can be checked with `python3 --version`.

This package depends on the NumPy and SciPy libraries. The examples additionally use Matplotlib for plotting. 

If you have not yet set up a python environment with NumPy and SciPy, this is easily accomplished on most platforms using [anaconda](https://docs.anaconda.com/anaconda/install/).


### About
ltkinetics was written by Zachary Mathe, doctoral student in the [DeBeer Group](https://cec.mpg.de/1/research/1087/prof-dr-serena-debeer/). 

This is a work in progress and feedback is welcome. Though I verified the implimentation by duplicating figures from Wilson's thesis (see example 4), users are encouraged to read through main.py and the equations therein. In the future, I may add functionality for reoptimizing kinetic constants and/or starting conditions against experimental data.

In the future, this package will be associated with either a Zenodo DOI or a peer-reviewed publication. If you use this package in your own work, please check back here for something to cite.


### References
1. [Wilson's Thesis, advised by Watt](https://scholarsarchive.byu.edu/etd/516/)
