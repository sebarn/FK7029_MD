The aim of this project is to produce stable computer simulations of a Lennard-Jones Argon liquid by using Molecular Dynamics and Monte Carlo techniques. Once the simulation code runs stable, measurements of respective energies and the heat capacity are performed. Moreover the radial distribution function is computed. The foundation of this project are the somewhat vague lab instructions found in \cite{LabInstructions}.

A personal side-challenge of this project is to write my first simulation using Python and trying to refute the common misconception that Python may not be used to write physics simulations. To run this code you need to install the  following:

* Python 2.7.3
* SciPy	(for computations)
* Numpy (for computations)
* vPython (for 3D visualization)
* matplotlib (for plotting in data analysis part)	

This code has been tested on Linux Ubuntu 12.04 LTE, 32-bit. The reason for installing Python 2.7.3 (instead of the latest version) is that the package SciPy.Weave currently only runs under Python 2.7. All of the above can be installed via the Ubuntu software center. In this code I try to find the right balance between readability, taking advantage of Python features and performance.
