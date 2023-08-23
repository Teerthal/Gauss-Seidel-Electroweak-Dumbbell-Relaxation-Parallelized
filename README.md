# Gauss-Seidel-Electroweak-Dumbbell-Relaxation-Parallelized
Parallelized Fortran routines to compute the electroweak dumbbell configuration using initial guess profiles and Gauss-Seidel numerical relaxation.

-----------------------------------------------------------------------------------------------------------------------------------------------------

# Description:
---------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------

These routines can fully reproduce the results in https://arxiv.org/abs/2302.04886

The purpose of this numerical suite is to compute the field configuration of the electroweak dumbbell described in the aforementioned article. 

This was accomplished by using initial guess profiles and using numerical relaxation, subject to the static electroweak equations of motion.

-----------------------------------------------------------------------------------------------------------------------------------------------------

## Code Description and limitations:
---------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------

The code computes the 3-dimensional relaxed field configuration of the Higgs scalar field, and the W and Y gauge fields in the standard bosonic electroweak 
theory. The numerical relaxation is made by iterating over the spatial points, with each update determined by the equations of motion stated in the reference above.

## Modules:
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

### LIMITATIONS:

-Needs a minimum of 2 processes in for each spatial dimension

-The interval between iterations for computing energies should be larger than N^3, where N is the number of processes in a single spatial dimension

-The number of processes should be cubes of integers. This can be remedied but has not been done yet.

## OUTPUT DATA:
-----------

The primary output is the total energy of the configuration at stated intervals until the simulation ends

Commented lines in main, energy and electromagnetics can enable output of |\Phi|, energies and magnetic fields as a function of xz in the xz plane(y=0)

The python scripts can read the stated format in these commented lines and produce xz-contour plots.

# Parallelization:
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

-----------------------------------------------------------------------------------------------------------------------------------------------------
## PARALLELIZATION SCHEME
We describe the parallelization scheme here:

- Divide the physical lattice into sublattices with appropriate ghost zones and assign them to unique parallel processes.
- (Step 1)Compute the first iteration in the first process and communicate with neighboring processes. Then immediately move on to the second iteration.
- (Step 2)The neighboring processes receive the communication and compute the first iteration, immediately followed by them performing the relevant communication, and moving to the second iteration.

-This would be repeated until all the processes are busy. If there are N sublattices(processes) in a single direction, then by the time the Nth process performs the 1st iteration, the 1st process would be performing the Nth iteration.

The advantage of this approach is that once all the processes have been saturated, there will be no free or unused processes until the time when global integration or data collection needs to be performed.

###Example
It is best to follow a concrete example to understand our parallelization implementation
We will use a 2-Dimensional lattice that is decomposed into 8^2 sublattices or subdomains, each on some process n.
Consider the process topology to be determined by coordinates (I,J) as illustrated here:
##Process Topology
![processortopology](https://github.com/Teerthal/Gauss-Seidel-Electroweak-Dumbbell-Relaxation-Parallelized/assets/95438989/e997a4d4-eff1-4e44-99a3-f581fa8a02bc)

##Step1: Conduct first iteration (it:0) over the lattice points in the first process (0,0) and send the updated field values to the ghost zone or halo in the processes (0,1) and (1,1) as illustrated below.
![it1](https://github.com/Teerthal/Gauss-Seidel-Electroweak-Dumbbell-Relaxation-Parallelized/assets/95438989/85579d8f-40dd-4b11-811e-74e04e1d99ca)

##Step2: 
- Process 1 (0,0) moves to compute iteration (it:1).
- Processes (1,0) and (0,1) receive communication from Process (0,0) and compute iteration (it:0)
- Processes (1,0) and (0,1) send communication to its neighboring process (Indicated by Blue arrows) [!Important Note: The horizontal arrows represent sending/receiving information in the horizontal direction and the vertical arrows represent the same but for the vertical direction. Therefore, there is no issue in different processes sending and receiving information to a common process, as long as it is done in different direction. Ex: Process(0,1) sending and receiving vertical halo information with (1,1) while Process(1,0) sends and receives horizontal halo information with Process(1,1).]
  - Process (0,1) -> (0,0),(1,1) & (0,2)
  - Process (1,0) -> (0,0),(1,1) & (2,0)

![it2](https://github.com/Teerthal/Gauss-Seidel-Electroweak-Dumbbell-Relaxation-Parallelized/assets/95438989/4bd7e135-a925-4a16-8bd1-a68edb031c80)

##Step 3: Same as the previous step but 1 iteration further along
![it3](https://github.com/Teerthal/Gauss-Seidel-Electroweak-Dumbbell-Relaxation-Parallelized/assets/95438989/bf437acd-1527-4ce8-87b0-9c69b6964f44)
.
.
.
.
##Step 8: The moment the first process (0,0) is computing the 8th iteration (it:7), the 8th process along the horizontal direction (0,7) will be conducting the first iteration (it:0).
![itn](https://github.com/Teerthal/Gauss-Seidel-Electroweak-Dumbbell-Relaxation-Parallelized/assets/95438989/dbe10bb7-cb91-4835-a514-3ef1be6609d4)


###Key Ingredient:

The key ingredient that makes this work in the code is that every communication conducted by the process is passed through MPI Send with a message (msg:it) marked by the iteration. This allows the parallel processes to function asynchronously and no process would need to wait after sending their share of the relevant data, and continue computing the next iteration. 

This, of course, provides maximal occupation of the processes. However, it comes with the caveat that, any MPI Barrier processes be only conducted after intervals between iterations that are greater than N^3 for 3 dimensions. This is because upon executing a barrier, say for example, to compute the total energy in the global domain at some iteration, all processes would need to reach the same iteration. Therefore, if the interval is less than N^3, this would imply that the entire process grid is never fully saturated. This is still advantageous over serial computing but the efficiency drops and would be the same (or even worse if taking into MPI Send/Receive overhead into account) if the interval = 1.
