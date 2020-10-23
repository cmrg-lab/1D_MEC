#Copyright © 2020 Viviane Timmermann, Simula Research Laboratory, Norway

## Timmermann V, McCulloch AD (2020). Mechano-electric coupling and arrhythmogenic
#  current generation in a computational model of coupled myocytes.
# Frontiers in Physiology, doi:

---

### What is this repository for? ###
The repository contains all data sets and computational code for isolated and
coupled cell simulations in the folders 'SingleCell' and 'MultiCell', respectively.
In the folder 'MultiCell' the data sets and code are divided in two different
folders for 3 and 6 coupled cell simulations ('3CoupledCells' and '6CoupledCells',
respectively).
Each of the folders 'SingleCell', '3CoupledCells', and '6CoupledCells' includes
sub-folders for 0, 1, and 3 coupled fibroblasts ('Fib0', 'Fib1', and 'Fib3',
respectively). These folders contain the conditions for:

* Sarcomere length random distribution (SarcRest*.txt)
* Fibroblast localisation and conductance random distribution (G_fib*.txt)
* Initial conditions of the ODE system (Init*.txt)

The code for the the different stretch patterns H, I, P, and S correspond to
the c++ models MultiMyocyteHealthy.cpp, MultiMyocyteIsometric.cpp,
MultiMyocytePlateau.cpp, and MultiMyocyteSinusoidal.cpp, respectively.

The variable 'cells' in l. 11 defines the number of coupled cells
(for isolated cell simulations cells = 1)

The variable 'iCount' in l. 1331, 1340, and 1357 defines the myocyte conditions
for which the simulation is run (Data Supplement of the manuscript in Tables 1,2,3)


### How do I run the simulation? ###

'MyocyteSim' has to be substituted by 'MultiMyocyteHealthy', 'MultiMyocyteIsometric',
'MultiMyocytePlateau', or 'MultiMyocyteSinusoidal'.

* g++ MyocyteSim.cpp -o make
* ./make
  or ./make >> datafile.txt (to save the output in 'datafile.txt')


### Who do I talk to? ###

* Viviane Timmermann, viviane@simula.no
* Andrew D. McCulloch, amcculloch@ucsd.edu
