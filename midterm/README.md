# CEE 576 Nonlinear Finite Element - Midterm Exam

## File descriptions

|Files   | Provided    | Modified | Category | Type | Notes |
|---------|-------------|-------------|-------------|------------|------------|
|`plot2DLagransurf.m`       | &#9989;   | &#10060; | Graphics | Function | |
|`plot2DLagransurfcont.m`   | &#9989;   | &#10060; | Graphics | Function | |
|`plotModel.m`              | &#9989;   | &#10060; | Graphics | Function | |
|`plotModelCont.m`          | &#9989;   | &#10060; | Graphics | Function | |
|`LagrSurfacePoint.m`       | &#9989;   | &#10060; | Graphics | Function | |
|`LagrSurfacePointCont.m`   | &#9989;   | &#10060; | Graphics | Function | |
|`LocToGlobDOF.m`           | &#9989;   | &#10060; | FEM helper function | Function | Map equation indices to node numbering |
|`shpl_2d.m`                | &#9989;   | &#10060; | FEM helper function | Function | Return 2D shape functions |
|`shpg_2d.m`                | &#9989;   | &#10060; | FEM helper function | Function | Return 2D shape gradients |
|`intpntq.m`                | &#9989;   | &#10060; | FEM helper function | Function | Return quadralateral quadrature |
|`intpntt.m`                | &#9989;   | &#10060; | FEM helper function | Function |  Return triangular quadrature |
|`GaussPoints.m`            | &#9989;   | &#10060; | FEM helper function | Function | Return Gauss-Legendre points |
|`Elas2d_Elem.m`            | &#9989;   | &#10060; | Element type def | Function | 2D Linear elastic element|
|`L_Elem2_2d.m`             | &#10060;  | &#9989;  | Element type def | Function | 2D Small strain nonlinear element |
|`AssemStiff.m`             | &#9989;   | &#10060; | Assembly | Script | 
|`AssemStiffForce.m`        | &#9989;   | &#10060; | Assembly | Script |
|`FormFE.m`                 | &#9989;   | &#9989;  | Assembly | Script | Assemble the global stiffness matrix
|`triangtwo.m`              | &#9989;   | &#10060; | Input | Script | Triangular element compression
|`input_quad_stretch.m`     | &#10060;  | &#9989;  | Input | Script | Quadrilateral element stretch (4a)
|`input_quad_shear.m`       | &#10060;  | &#9989;  | Input | Script | Quadrilateral element shear (4b)
|`assign_bc_load_data.m`    | &#9989;   | &#10060; | FEM core routine | Script | 
|`CompStrainStress.m`       | &#9989;   | &#10060; | FEM core routine | Script | Compute global strain and stress
|`CompStrainStress_Elem.m`  | &#9989;   | &#10060; | FEM core routine | Function | Compute elemental strain and stress
|`FintStressStrain.m`       | &#10060;  | &#9989;  | FEM core routine | Script | Compute nonlinear strain, stress, and F_int
|`SolveFE.m`                | &#9989;   | &#9989;  | Solve | Script | 
|`FEA_Program.m`            | &#9989;   | &#9989;  | Driver | Script |
|`load_driver.m`            | &#10060;  | &#9989;  | Driver | Script | Main driver file for Problem 4

## Instruction for running Problem 4 source code
To reproduce the figures for Problem 4, simply run in MATLAB `load_driver.m`. 
The boolean variable `q4a` toggles between problem 4(a) and 4(b). If `q4a` is zero, case 4(b) (shear) is chosen, otherwise 4(a) (stretch) is chosen. 
The parameters are set in the two input files `input_quad_stretch.m` and `input_quad_shear.m`. 