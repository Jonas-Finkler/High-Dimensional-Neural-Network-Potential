# High-Dimensional-Neural-Network-Potential

Fortran implementation of a high dimensional neural network potential.

This is a simple implementation of the prediction step for a high dimensional neural network potential (HDNNP)[1].
Training of the HDNNP is not implemented. 
Subroutines are included to parse the parameters fitted by RuNNer [2].
However only the most basic features are implemented and only the G2 and G4 symmetryfunctions [3] can be calculated. 
The code is designed in a modular way, such that the implementation of additional features should be simple.


[1] Behler, Jörg, and Michele Parrinello. "Generalized neural-network representation of high-dimensional potential-energy surfaces." Physical review letters 98.14 (2007): 146401.  
[2] [https://www.uni-goettingen.de/de/560580.html](https://www.uni-goettingen.de/de/560580.html)  
[3] Behler, Jörg. "Atom-centered symmetry functions for constructing high-dimensional neural network potentials." The Journal of chemical physics 134.7 (2011): 074106.  

