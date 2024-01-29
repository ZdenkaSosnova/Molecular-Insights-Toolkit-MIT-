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

## Methods of the Class
The `Bond_lenght_Analyzer` class provides the following methods:

1. **`__init__(file_xyz, dimension=2, minimal_distance=1.35, maximal_distance=2)`**: Initializes the analyzer with molecule coordinates from a specified file (`file_xyz`). Users can customize the minimum and maximum bond lengths (`minimal_distance` and `maximal_distance`) as needed. The dimensionality of the molecule can be specified (2D or 3D) using the `dimension` parameter.

2. **`graph_2d()`**: Generates a 2D graph depicting bond lengths in the molecule defined by coordinates in ".xyz" format. Different bond lengths are represented with different colors, and each bond is labeled with its length.

3. **`projection_y_z_axis()`**: Generates 3D projections of bond lengths analysis viewed from the Y-axis and Z-axis.

4. **`graf_3d()`**: Generates a 3D graph depicting bond lengths in the molecule.


## Exciting Possibilities
The `Bond_lenght_Analyzer` class opens up exciting possibilities for studying molecular structures and properties. By visualizing and analyzing bond lengths, researchers and students can gain deeper insights into the intricate world of chemistry and molecular biology. Whether exploring fundamental principles or tackling complex research problems, the `Bond_lenght_Analyzer` class is a versatile tool for unraveling the mysteries of molecular systems.
