# Huckel Method Implementation (huckle_method.py)

## Overview

The `huckle_method.py` module provides implementations of the classical and extended versions of the Huckel method for analyzing the electronic structure of planar, pi-conjugated molecules. This module aims to compute and visualize eigenvalues and eigenvectors, offering users valuable insights into molecular orbitals and electronic properties.

## Introduction to the Huckel Method

The Huckel method, named after Erich Hückel, is a powerful computational approach used to study the electronic structure of conjugated organic molecules. It provides approximate solutions to the Schrödinger equation for pi-electron systems, focusing on planar molecules with delocalized pi-electrons. The method is based on the tight-binding approximation and involves parameters such as the on-site energy (alfa) and the interaction between neighboring atoms (beta). By solving the Huckel matrix, which represents the system's Hamiltonian, the method yields eigenvalues and eigenvectors corresponding to molecular orbitals and their energies.


### Installation and Necessary Libraries

Ensure you have Python 3.x installed.
The `huckle_method.py` module relies on the following libraries:
- **matplotlib**: Used for plotting graphs and visualizing molecular orbitals.
- **numpy**: Essential for numerical calculations and matrix operations.
- **scipy**: Utilized for advanced scientific computing tasks such as solving eigenvalue problems.

## Methods of the Class

### Constructor (`__init__`)

- **Parameters**:
  - `file`: File in ".xyz" format specifying the coordinates of the molecule.
  - `alfa`: On-site energy of the atom.
  - `beta`: Parameter determining the interaction between individual atoms.
  - `extended_huckel`: Boolean variable indicating whether to use the extended Huckel method.
  - `number_of_states`: Number of states around the Fermi energy to be represented.
  - `minimal_distance`: Minimum distance between individual atoms.
  - `maximal_distance`: Maximum distance between individual atoms.

### Additional Methods

- `energy_graph()`: Generates a plot of energies (eigenvalues) around the Fermi energy.
- `orbital_graph(orbital, state)`: Creates a graph for the selected orbital.
- `huckel_orbitals()`: Plots selected molecular orbitals around the Fermi energy.
- `return_gap_value()`: Calculates the energy difference between the highest occupied and lowest unoccupied orbital.
- `bond_charge_matrix_txt()`: Saves the bond charge matrix to a text file.
- `graph_bond_charge()`: Generates a visualization of the bond charge matrix.

## Time Characteristics Analysis

An additional analysis was conducted to assess the computational time required for molecules ranging from a number of atoms to 2000. All of these times were found to be lower than a linear function 0.004x + 0.002, where x represents the number of atoms. It's important to note that this analysis is dependent on the user's computer and was conducted for a specific setting including plotting energy graphs and plotting 4 molecular orbitals. The purpose of this analysis is to provide users with an idea of the computational time range, given the square matrices' shape defined by the number of atoms.
![Graph of "time analysis"](https://previews.dropbox.com/p/thumb/ACJWoYplvAcxCgwO2aGuIYqcvumkjPwhSh0l58x84UikyANrSOOctmEnszeBRXiUjuojm9mbIHe1n_U-p2zlh1qm3vl8sc8m-erGJfxK_XckPvXZBSoHea_B50e4knOFSDt1pfKq3RMwMO3Ty4GY3fmoi8K2NPHFJBjodV1XlgAcaVfgSaB74VCHnwxz17V_As3_UItPNo6bmM07Tmmp9uld8O8PX54HTmpNOntLr8aGvbKqys68ISVhA6v3F9N7DgUWu2aH_rHksH7AtPtpQA5qm3VwJZQng1ACQmNlGI9rZtSle6ZaNKsuIHCMxTLyhVroODcgliiAV66jgYitQLXl/p.png)




