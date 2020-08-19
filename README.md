## Ising Decimation Reversion with CNN

Metropolis Simulation of the 1D Ising model. Reverting a renormalization by decimation with a trained convolutional neural network (CNN). 

## Background

The Ising model is a nearest neighbor interaction approximation of a magnet with each spin being arranged in a binary state on a lattice. It serves as a very important toy model in computational physics and there are many noticeable problems with it which also apply to more complex models. One example is the critical slowing down in computational simulations of the Ising model at the critical temperature of the phase transition. It is a problem which is originated from a divergent correlation length of the interactions. This is especially crucial when studying larger system sizes, since it would take very long to capture all important states in a more traditional computer simulation, like a Monte Carlo (MC) simulation with a Metropolis algorithm.

When using a renormalization procedure on a Ising spin configuration to decrease the system length to a half, the idea is then to train a neural network with supervised learning to revert the renormalization to its original system size. It is not required for the super resolution to be an exact reversion to the problem, since only the probability distribution is what plays a role when calculating a statistical average of an observable. It would then be possible to run a computer simulation with a MC method of the Ising model at small system sizes and rescale their properties to a much larger size. The idea comes from a paper from [Efthymiou, Beach and Melko, which was published on the 31st of January 2019](https://arxiv.org/abs/1810.02372).

This work is a groundwork for an enlarging method. The focus is on reconstructing a decimation procedure to the original configuration distribution.

## C++ Metropolis Algorithm

The Metropolis simulation is used to create a dataset of multiple spin configurations. This is executed with the C++ code in 'ising.cpp'. C++ is chosen for efficiency reasons. The spin configurations are also decimated and enlarged by the majority rules in this code. The spin configurations are saved as a .csv file. Three observables (absolute magnetization, energy and spin-spin correlation function) are also calculated here.

## Python

The CNN reconstruction is learned with the code in 'learning.ipynb'. A CNN model is prepared and trained to convert the decimation procedure to match the original spin configuration distribution. Then it is shown, that this reconstruction does match the correct observables.

## Contribution

The code is written by Jan Zimbelmann
