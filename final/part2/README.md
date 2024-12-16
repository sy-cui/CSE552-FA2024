# CEE 576 Final Exam Part II

## Outline of the `final/part1` directory
* Files at the top level are main driver and execution files for the FEM program. 
* `visual/`: scripts and routines that generate visualizations.
* `input/`: all input files are located in this directory.
* `utility/`: generic helper functions for FEM.
* `elasticity`: routines and definitions pertaining specific constitutive relations. Files related to this exam are all located with the `elasticity/hyperelasticity/` directory. 
* `dynamics/`: helper functions for dynamical analysis, such as routines for computing predictors and correctors for the HHT-method.

## Important new additions
* All files in `dynamics/`. 
* `utility/FormMass.m`: script that assembles the mass matrix after calling `dynamics/mass_elem_2d.m`
* `Fint.m`: a master file for computing and assembling the internal force vector. This file is structured very similarly to `FormFE.m` except that the latter forms the consistent tangent matrix. 
* New structure for `elasticity`: every subdirectory (e.g. `elasticity/linear/`) nows contain a `fint_elem_...` and a `tang_elem_...` function, which computes the element internal force and consistent tangent, respectively. These routines are called within `Fint.m` and `FormFE.m`.
* `elasticity/.../ss_..._2d.m`: stress and strain computation for every constitutive relation. Typically only called during post-processing should these properties be desired. 

## Input files 
The two input files related to this exam is located at `input/input_stretch.m` and `input/input_shear.m`. 

## Executing the program
Case (b) is executed by `dynamic_stretch.m`, which runs `input/input_stretch.m` internally. Case (c) is executed by `dynamic_shear.m`, which runs `input/input_shear.m` internally. 

Please make sure that `header.m` is executed before any other scripts. Otherwise MATLAB does not load the subdirectories into its path.