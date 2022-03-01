# High-Dimensional-Neural-Network-Potential

Fortran implementation of a high dimensional neural network potential.

This is a simple implementation of the prediction step for a high dimensional neural network potential (HDNNP)[1].
Training of the HDNNP is not implemented. 
Subroutines are included to parse the parameters fitted by RuNNer [2].
However only the most basic features are implemented and only the G2 and G4 symmetryfunctions [3] can be calculated.
Electrostatics (3G- and 4G-HDNNP) are also not implemented. However, the implementation should not be too complicated using my [Ewald summation code](https://github.com/Jonas-Finkler/ewald-summation).
A standalone version of the neural network implementation including a training example using the Kalman filter can be found [here](https://github.com/Jonas-Finkler/fortran-NeuralNetwork). 

The code is designed in a modular way, such that the implementation of additional features should be simple.

A simple example running a molecular dynamics simulation on an example HDNNP for a methylammonium molecule is included in the file ```src/potentials/hdnnpExample.f90```.
The parameters for the potential can be found in the ```example``` directory.

Atomic units are used for all variables in the code. The calculation of the symmetryfunctions is parallelized using openmp.
In most cases the evaluation of the symmetryfunction is the most expensive part. 
However, for large systems the construction of the neighbourlists might become more expensive since the current implementation 
loops over all possible pairs of atoms. 
In that case it might be advantageous to implement an algorithm that bins the atoms into a grid.
Since the evaluation of the neural networks is usually rather inexpensive, calls to BLAS are not yet implemented. 

A CMake script to compile the code is included. It can be used as follows:
```bash
mkdir build
cd build
cmake .. -DINTEL=OFF -DDEBUG=OFF -DOPENMP=ON #-DINTEL=OFF uses gfortran =ON uses ifort
make 
cd ../example
../build/hdnnp.x # run the example program
```

[1] Behler, Jörg, and Michele Parrinello. "Generalized neural-network representation of high-dimensional potential-energy surfaces." Physical review letters 98.14 (2007): 146401.  
[2] [https://www.uni-goettingen.de/de/560580.html](https://www.uni-goettingen.de/de/560580.html)  
[3] Behler, Jörg. "Atom-centered symmetry functions for constructing high-dimensional neural network potentials." The Journal of chemical physics 134.7 (2011): 074106.  

