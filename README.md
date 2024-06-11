# Molecular-Dynamics-Simulation-of-Particles-in-a-Box
A C++ molecular dynamics simulation that models particle behaviour in a 3D box using the Lennard-Jones potential and the Metropolis algorithm.

This project implements a basic molecular dynamics simulation to model the behavior of particles within a three-dimensional box. The simulation initializes particle positions on a lattice and uses the Lennard-Jones potential to compute the energy of the system. Over a series of steps, particles are randomly displaced, and the Metropolis algorithm is employed to determine whether the new positions are accepted based on the change in energy. The simulation outputs the positions of the particles at regular intervals, which can be visualized using tools that support XYZ format.


Features:

-Initializes particles on a lattice within a cubic box.

-Computes system energy using the Lennard-Jones potential with a cutoff.

-Randomly displaces particles and uses the Metropolis algorithm for acceptance.

-Outputs particle positions at regular intervals in XYZ format.


Requirements:
C++11 or later


Usage:

-Compile the code using a C++ compiler, e.g., g++ -o simulation simulation.cpp.

-Run the executable, ./simulation.

-The particle positions are saved in pos_eps10e-1.xyz file, which can be visualized using compatible visualization tools.


Detailed Description of Key Functions:

-initializeParticles: Sets up initial particle positions on a lattice.

-computeEnergy: Calculates the total potential energy of the system using the Lennard-Jones potential.

-updatePositions: Randomly displaces a selected particle within the box.

-initialParticlesPosition: Outputs the initial positions of particles.

-write: Writes the current positions of particles to a file at specified intervals.
