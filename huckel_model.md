# Huckel Method Implementation (huckle_method.py)

## Overview

The `huckle_method.py` module provides implementations of the classical and extended versions of the Huckel method for analyzing the electronic structure of planar, pi-conjugated molecules. This module aims to compute and visualize eigenvalues and eigenvectors, offering users valuable insights into molecular orbitals and electronic properties.

## Introduction to the [Hückel Method](https://daniloroccatano.blog/2018/05/23/the-simple-huckel-method/)

The Huckel method, named after Erich Hückel, is a powerful computational approach used to study the electronic structure of conjugated organic molecules. It provides approximate solutions to the Schrödinger equation for pi-electron systems, focusing on planar molecules with delocalized pi-electrons. The method is based on the tight-binding approximation and involves parameters such as the on-site energy (alfa) and the interaction between neighboring atoms (beta). By solving the Huckel matrix, which represents the system's Hamiltonian, the method yields eigenvalues and eigenvectors corresponding to molecular orbitals and their energies.


### Installation and Necessary Libraries

Before using the `Huckel model` class, ensure you have installed the following Python libraries:
- `utils.py`: A utility library for handling molecule coordinates and distance calculations.
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
- `bond_charge()`:bond_charge(): This function calculates the bond charge, which measures the strength of a pi-bond between any two atoms (i,j) in two molecules. For further analysis, only the strength of bonds between nearest neighbor atoms is considered. You can learn more about bond charge [here](https://www.chm.bris.ac.uk/pt/ajm/html/L4_p2.htm).
- `bond_charge_matrix_txt()`: Saves the matrix to a text file. 
- `graph_bond_charge()`: Generates a visualization of the bond charge matrix.

## Time Characteristics Analysis

An additional analysis was conducted to assess the computational time required for molecules ranging from a number of atoms to 2000. All of these times were found to be lower than a linear function 0.004x + 0.002, where x represents the number of atoms. It's important to note that this analysis is dependent on the user's computer and was conducted for a specific setting including plotting energy graphs and plotting 4 molecular orbitals. The purpose of this analysis is to provide users with an idea of the computational time range, given the square matrices' shape defined by the number of atoms.
<p align="center">
  <img src="https://ucbcee984ddbed20fe99af091fc7.previews.dropboxusercontent.com/p/thumb/ACLtGvnpA278zUZSA3wh9lcxVaGmM1RzxW27cWErVPeL0MjGSuM7Fc_4Nu0VbIKLR79mdrIDZQ6NSeg2F_u34yovyX6I5igSyrRUBrhU3FgJAE0tNs-zz5JBxftgn1AvS9ukqyR1wTQ4Vej-RWA3WbauyanWjiAkB4imdwByqpDWw6eIJ1n-7pMGpANlzQg86frK8zYrK-_NarxuuukyAu-zIbh4sdRNP3DgjXnPZHmkEBpkOKHSpEFRcWZ8G0pZri2B69WYM2bmA1IiYWxBRxDD0lXwjfTvx3_KbopvQeSPZ31Jo20QIfBXMSGrGcXq1vbyWHMOD6Nd3qpwpWmFBpodC8T_gaB0V_VSxkq4QfnOdV4vb7trcuG0TdXYMd2Gd1Y/p.png" alt="Graph Title">
</p>






