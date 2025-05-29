# MoS2Nanopores
Cataloging isomers of nanopores (extended vacancy defects) in 2D lattices of TMDs

> Bhowmik, S.; Govind Rajan, A. Size and Chemical Environment Influence Structural Order in Defective MoS<sub>2</sub> from Irregular to Triangular Nanopores

## Contents

**All codes are written in MATLAB 2024a. The following are included programs:**

* `generate_isomers_SiEtch`: This calls the `kmc_isomers_SiEtch` repeatedly to generate a pre-determined number of isomers. In so doing, the program requires:
    * `kmc_isomers_SiEtch`:This program carries out a kinetic Monte Carlo (KMC) simulation to stochastically generate one nanopore of a given size in the graphene lattice, and saves the XYZ file and the directed adjacency matrix of each generated nanopore.

    * `barrier_calc`: The `kmc_isomers_SiEtch` program calls this to calculate the barrier of each atom etching event according the atomic fingerprints.

* `analyze_isomers_directed_isomorphism`: This program analyzes all the nanopore isomers generated using the `generate_isomers_SiEtch` code to weed the duplicate isomers out of the list, and output the unique isomers. To do so, the program requires:

  * `add_weighing_nodes_in_between`: This program modifies the antimolecule adjacency matrix to add fictitious nodes to enable proper detection of C-C bonds with various orientations 

  * `compare_pores_with_weights`: This program compares two nanopore isomers, after fictitious nodes have been added into their adjacency matrices using `add_weighing_nodes_in_between`, and outputs whether the two nanopores are identical or not. 
