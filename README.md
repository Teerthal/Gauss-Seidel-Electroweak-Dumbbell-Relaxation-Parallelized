# Gauss-Seidel-Electroweak-Dumbbell-Relaxation-Parallelized
Parallelized Fortran routines to compute the electroweak dumbbell configuration using initial guess profiles and Gauss-Seidel numerical relaxation.

-----------------------------------------------------------------------------------------------------------------------------------------------------

Description:
---------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------

These routines can fully reproduce the results in https://arxiv.org/abs/2302.04886

The purpose of this numerical suite is to compute the field configuration of the electroweak dumbbell described in the aforementioned article. 

This was accomplished by using initial guess profiles and using numerical relaxation, subject to the static electroweak equations of motion.

-----------------------------------------------------------------------------------------------------------------------------------------------------

Code Description and limitations:
---------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------

The code computes the 3-dimensional relaxed field configuration of the Higgs scalar field, and the W and Y gauge fields in the standard bosonic electroweak 
theory. The numerical relaxation is made by iterating over the spatial points, with each update determined by the equations of motion stated in the reference above.

Modules:
-------

boundaryIndices,processCoord: These subroutines determine the process topology for parallelization

ghostrecvl,ghostrecvu,ghostsend: These subroutines enable communication between neighboring sublattices.

main: Main file that calls all the subouines, instantiates the initial guess configuration, carries out Gauss-Seidel relaxation, and collects energies and magnetic field data at specified intervals.

derivatives4thOrder: Numerical differentiation scheme that uses 6th order expressions in the lattice interior and 2nd order derivatives in the 3 boundary layers at the edge of the lattice.

energy: Computes and integrates various energy components.

fluxesNumRel: Computes the fluxes determined by the equation of motion and used by the Euler revolver

evolveeuler: Computes the updated field values using the fluxes and carries out communication with neighboring lattices

evolveloop: Called by evolveeuler to iterate over the spatial points within the relevant sublattice and impose boundary conditions

covariantDerivs: Computes the covariant derivatives of the Higgs field

fieldStrengths: Computes the gauge field strengths

printenergy: Collects energies from all the domains on the 0th process, adds them and then prints them

icMMbar,icgauge: Instantiates the initial guess field configuration and the Higgs configuration that is held fixed throughout the relaxation.

LIMITATIONS:

-Needs a minimum of 2 processes in for each spatial dimension

-The interval between iterations for computing energies should be larger than N^3, where N is the number of processes in a single spatial dimension

-The number of processes should be cubes of integers. This can be remedied but has not been done yet.

OUTPUT DATA:
-----------

The primary output is the total energy of the configuration at stated intervals until the simulation ends

Commented lines in main, energy and electromagnetics can enable output of |\Phi|, energies and magnetic fields as a function of xz in the xz plane(y=0)

The python scripts can read the stated format in these commented lines and produce xz-contour plots.

Parallelization:
----------------
For the purpose of lattice simulation, it is common to divide the physical domain into subdomains handled by the different processes.
The division of the physical domains has to be accompanied by ghost sites or a halo around the physical domains which has information about the neighboring 
domains. This is required for various reasons, the most common one, and the one relevant in this code, is to compute numerical derivatives.

The most common approach in studies of evolution would be:

-Compute a time step across each process

-Collect any data across all the domains required for computing global variables such as the total energy in the simulation box

-Update the halo around each domain with information from neighboring domains

-Compute the next step and repeat until some criteria have been met.

This is not the most efficient method for the Gauss-Seidel relaxation method.

For Gauss-Seidel relaxation, as one iterates over the lattice sites, the updated information from the previous lattice site is immediately available
for computing the updated field values. This is advantageous over Jacobi method since it has been shown to have faster convergence.

This is a challenge since a single domain has to be updated and its halo information has to be sent to the neighboring domains for them to begin computation,
and cannot be achieved by the simple domain decomposition scheme described above.

An efficient parallelization scheme has been implemented here, built upon the framework by Aayuush Saurabh and used for the works:

https://arxiv.org/abs/1705.03091

https://arxiv.org/abs/1904.02257

We describe the parallelization scheme here:

-Divide the physical lattice into sublattices with appropriate ghost zones and assign them to unique parallel processes.

-Compute the first iteration in the first process and communicate with neighboring processes. Then immediately move on to the second iteration.

-The neighboring processes receive the communication and compute the first iteration, immediately followed by them performing the relevant communication, and moving to the second iteration.

-This would be repeated until all the processes are busy. If there are N sublattices(processes) in a single direction, then by the time the Nth process performs the 1st iteration, the 1st process would be performing the Nth iteration.

The advantage of this approach is that once all the processes have been saturated, there will be no free or unused processes until the time when global integration or data collection needs to be performed.
