% Arbitrary data for assistance in defining the mesh
L = 1;  % [m]
H = 1;  % [m]

% Mesh Nodal Coordinates
NodeTable = [0 0
             L 0
             0 H
             L H];
numnp = length(NodeTable);

% Current nodal coordinate. To be updated
NodeCurr = NodeTable;

% Mesh Element Connectivities
ix = [1 2 4 3 1];
nen = 4;
numel = 1;

% Mesh Boundary Conditions and Loads
BCLIndex = [6 0]';
NodeBC = [1 1 0
          1 2 0
          2 1 0
          2 2 0
          3 2 0.2
          4 2 0.2];
NodeLoad = [3 2 0.2
            4 2 0.2];

% Mesh Material Properties
young = 10;    % [MPa]
pois = 0.25;
lame1 = pois * young / (1 - 2*pois) / (1 + pois);
lame2 = young / 2 / (1 + pois);
thick = 1;
density = 1;

MateT = [lame1 lame2 thick density];

% Tolerance
tol = 1e-8;
max_iters = 50;

% Total / Updated Lagrangian
framework = 'UL';

% HHT time stepping
alpha = -1/3;
beta  = 0.25*(1-alpha)^2;
gamma = 0.5*(1-2*alpha);
dt = 1e-3;
tend = 10;
hht = [alpha beta gamma dt];

lumped = 0;
