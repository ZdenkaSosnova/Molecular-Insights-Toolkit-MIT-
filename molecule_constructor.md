# Molecule Constructor

## Overview
The `Molecule_constructor` class in the Molecular-Insight_Toolkit (MIT) Python library facilitates the construction of molecular structures by generating Cartesian coordinates based on specified parameters. This class is particularly useful for creating custom molecular configurations with precise control over bond lengths and arrangement. 

## Installation and Necessary Libraries
Before using the `Molecule_constructor` class, ensure you have installed the following Python libraries:
- `utils`: A utility library for calculating lengths and handling molecule coordinates.
- `matplotlib`: A plotting library for generating molecular graphs and visualizations.
- `numpy`: A numerical computing library for array manipulation and mathematical operations.

## Methods of the Class
The `Molecule_constructor` class provides the following methods:

1. **`__init__(length_inside, length_between, number_of_benzene_rings, repetition_count, file_name, **kwargs)`**: Initializes the molecule constructor with parameters such as bond lengths, benzene ring count, repetition count, and output file name. Additional elementary cell configurations can be specified using keyword arguments.

2. **`molecule_coordinates_to_xyz_file()`**: Writes the Cartesian coordinates of the generated molecule into a .xyz file format.

3. **`show_graph()`**: Displays a graph of the constructed molecule, visualizing bond lengths using different colors.

### Specific Molecules
The `Molecule_constructor` is optimized for creating molecules with the following characteristics:
- Only odd numbers of benzene rings per unit cell are supported.
- Any value of inside or between length can be specified.
- The generated `.xyz` file will only consist of carbon atoms (without hydrogen atoms).

## Example
Here's an example of how to use the `Molecule_constructor` class to generate a specific molecule:

```python
from Molecule_constructor import Molecule_constructor

# Create a Molecule_constructor instance with specified parameters
my_molecule = Molecule_constructor(length_inside = 1.23, length_between = 1.51, number_of_benzene_rings = 3,
                                   repetition_count = 12, file_name = "molecule", cell1=[1.25, 1.42, 5, [3,10]],
                                   cell2=[1.33, 1.21, 9, [8]], cell3=[1.55, 1.32, 1, [4]])

# Plot the molecule
my_molecule.show_graph()

```

## Exciting Possibilities
With the `Molecule_constructor` class, users can unleash their creativity in molecular design and exploration. Whether simulating novel organic compounds or studying crystal structures, this class offers a powerful tool for advancing research and education in chemistry and materials science. From simple organic molecules to intricate polymers, the possibilities are endless with the Molecular-Insight_Toolkit.
