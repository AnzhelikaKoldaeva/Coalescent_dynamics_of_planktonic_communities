# Coalescent_dynamics_of_planktonic_communities
The backward and forward models for modeling planktonic communities in chaotic advection. Examples to run the codes:

-------------------------------------------- Forward --------------------------------------------

#DIFFUSIVE
./Main -mutation_prob 0.1 -D_factor 0.1 -Ns 16000 -random_seed 1

#VORTEX
./Main -mutation_prob 0.1 -D_factor 0.1 -Ns 8192 -w_jet 1 -L_jet 0.25 -D_factor 0.1 -random_seed 1


-------------------------------------------- Backward --------------------------------------------

#DIFFUSIVE
./Main -mutation_prob 0.1 -D_factor 0.1 -Ns 16000 -random_seed 1

#VORTEX
./Main -mutation_prob 0.1 -D_factor 0.1 -Ns 8192 -w_jet 1 -L_jet 0.25 -D_factor 0.1 -random_seed 1


Parameters to modify:
-mutation_prob: Mutation probability
-D_factor: Diffusion factor
-Ns: Population size of the sample
-random_seed: random seed of the simulation (saved as file identifier too)

To consider advection flux:
Uncomment "#define ADVECTION" in define.h
