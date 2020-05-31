# Ising Model
Various numerical algorithms to solve the Ising model in 2 and 3D. These algorithms include the Metropolis algorithm, the Prokof'ev and Svistunov single worm algorithm. Also included is a more efficient extension to the Prokof'ev and Svistunov algorithm, which we called the double worm algorithm.

The code could definitely be written in a nicer, more modular way. However, each algorithm pulls from its corresponding class header file and from the "statistics.h" header file. The statistics header file contains basic statistics functions and some autocorrelation functions. 
