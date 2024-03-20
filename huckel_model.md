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

## Example
Here's an example of how to use the `Huckel model` class to generate a specific molecule:

```python
from molecule_constructor import Molecule_constructor
from huckel_model import Huckel_model

# Create file "molecule.txt"
my_molecule = Molecule_constructor(length_inside = 1.25, length_between = 1.55, number_of_benzene_rings = 5,
                                   repetition_count = 10, file_name = "molecule" )
my_molecule.molecule_coordinates_to_xyz_file()

# Calculate eigenvectors and eigenstates of the molecule
model = Huckel_model(file = "molecule.xyz", alfa = 0, beta = -2.8, extended_huckel = False, number_of_states = 4, minimal_distance = 1.20, maximal_distance = 1.60)
model.energy_graph()
model.huckel_orbitaly()


