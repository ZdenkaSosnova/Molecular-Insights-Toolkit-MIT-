# Bond Length Analyzer

## Overview
The `Bond_lenght_Analyzer` class in the Molecular-Insight_Toolkit (MIT) Python library is designed for analyzing bond lengths in a molecule based on relaxed coordinates. This class provides various methods for visualizing and interpreting bond lengths, offering valuable insights into molecular structures.

## Understanding Bond Lengths
The relaxation process of a molecule often results in changes to bond lengths between its constituent atoms. Analyzing bond lengths provides valuable insights into the stability and structural properties of molecules. Understanding bond lengths is crucial in various fields such as chemistry, materials science, and biochemistry. For example, in organic chemistry, precise knowledge of bond lengths aids in predicting molecular reactivity and behavior.

## Installation and Necessary Libraries
Before using the `Bond_lenght_Analyzer` class, ensure you have installed the following Python libraries:
- `utils.py`: A utility library for handling molecule coordinates and distance calculations.
- `matplotlib`: A plotting library for generating 2D and 3D visualizations.
- `numpy`: A numerical computing library for array manipulation and mathematical operations.


### Additional Methods

- `energy_graph()`: Generates a plot of energies (eigenvalues) around the Fermi energy.
- `orbital_graph(orbital, state)`: Creates a graph for the selected orbital.
- `huckel_orbitals()`: Plots selected molecular orbitals around the Fermi energy.
- `return_gap_value()`: Calculates the energy difference between the highest occupied and lowest unoccupied orbital.
- `bond_charge()`:bond_charge(): This function calculates the bond charge, which measures the strength of a pi-bond between any two atoms (i,j) in two molecules. For further analysis, only the strength of bonds between nearest neighbor atoms is considered. You can learn more about bond charge [here](https://www.chm.bris.ac.uk/pt/ajm/html/L4_p2.htm).
- `bond_charge_matrix_txt()`: Saves the matrix to a text file. 
- `graph_bond_charge()`: Generates a visualization of the bond charge matrix.

## Methods of the Class

### Constructor (`__init__`)

- **Parameters**:
  - `file_xyz`: File in ".xyz" format specifying the coordinates of the molecule.
  - `dimension`: dimensionality of the molecule can be specified (2D or 3D) using the `dimension` parameter
  - `minimal_distance`: Minimum distance betwween atoms, that should be visualized.
  - `maximal_distance`: Maximum distance between individual atoms, that should be visualized.
 
### Additional Methods

- **`graph_2d()`**: Generates a 2D graph depicting bond lengths in the molecule defined by coordinates in ".xyz" format. Different bond lengths are represented with different colors, and each bond is labeled with its length.

- **`projection_y_z_axis()`**: Generates 3D projections of bond lengths analysis viewed from the Y-axis and Z-axis.

- **`graf_3d()`**: Generates a 3D graph depicting bond lengths in the molecule represented with different colors.

