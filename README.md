# Monodisperse clustering external mixing (MCEM): A program for more realistic aggregation modeling of soot

The present code aims to generate fractal soot aggregates from clustering of primary particle populations. The aggregation phenomenon studied here is in the diffusion-limited cluster-cluster aggregation (DLCA) regime. The numerical algorithm used consists of solving transport equation for the particles in a Langeving Dynamics (LD) framewrok [1] and simulating their growth upon collision with each other. The primary particles considered are nascent soot in typical flame conditions having a polydisperse size distribution. The long-term goal for this project is to model the aggregation in a new way which assumes a hierarchical clustering suggested by External Mixing Hypothesis [2]. This will results in aggeragtes which are externally polydisperse (in terms of primary particle size) while having an internally monodisperse morphology.

Like many other LD-DLCA algorithms, MCEM follows the below steps to generate agglomerated soot:
1. Reading the aggregation input parameters,
2. Random initialization of the monomers' locations and velocities,
3. Moving the particles based on brownian motions and drag,
4. Imposing periodic conditions on the boundaries,
5. Checking for particle collisions over time,
6. Attaching the collided particles as larger clusters,
7. Maintianing the volume fraction by enlarging the main computational domain upon clusteraions,
8. Exporting and plotting the aggregation data.

To accomplish the above steps the following packages come with MCEM:
a. "PAR": Initializes different properties of particles.
b. "MOV": Solves the equation of motion for the particles and applies the proper boundary conditions.
c. "COL": Contains different tools to check collisions and cluster the particles.
d. "CPL": Assigns flow-particle couplings caused by interactions of soot with the background flaming environment.
In addition to these, the "input" and "output" folders, respectively, contain the initialization parameters and final results of the program.

References:
1. Suresh, V., & Gopalakrishnan, R. (2021). Tutorial: Langevin Dynamics methods for aerosol particle trajectory simulations and collision rate constant modeling. Journal of Aerosol Science, 155, 105746.
2. Olfert, J., & Rogak, S. (2019). Universal relations between soot effective density and primary particle size for common combustion sources. Aerosol Science and Technology, 53(5), 485-492.
