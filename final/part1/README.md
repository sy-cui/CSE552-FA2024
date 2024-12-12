# CEE 576 Final Exam Part I

## Outline of the `final/part1` directory
* Files at the top level are main driver and execution files for the FEM program, such as `FEA_Program.m` and `SolveFE.m`.
* `visual/`: scripts and routines that generate visualizations.
* `input/`: all input files are located in this directory.
* `utility/`: generic helper functions for FEM.
* `elasticity`: routines and definitions pertaining specific constitutive relations. Files related to this exam are all located with the `elasticity/hyperelasticity/` directory. 

## Input files 
The two input files related to this exam is located at `input/input_nh_stretch.m` and `input/input_nh_shear.m`. 

## Executing the program
Simply run in MATLAB `final1.m`, which load inputs from one of the two files. One can switch between the files by comment / uncomment line 5 and 6 in `final.m`. The load steps and increment can be modified in the same file, defined as array `P`. 

The uniaxial test program can be replicated by running `final1_test.m`. 

Make sure that `header.m` is executed before any other scripts. Otherwise MATLAB does not load the subdirectories into its path.

