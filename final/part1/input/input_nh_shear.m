% Arbitrary data for assistance in defining the mesh
L = 1;  % [m]
H = 1;  % [m]
% P = linspace(0, 1e7, 100); % [N]

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
BCLIndex = [4 2]';
NodeBC = [1 1 0
          1 2 0
          2 1 0
          2 2 0];
NodeLoad = [3 1 0
            4 1 0];

% Mesh Material Properties
young = 1e8;    % [MPa]
pois = .25;
lame1 = pois * young / (1 - 2*pois) / (1 + pois);
lame2 = young / 2 / (1 + pois);
thick = 1;

MateT = [lame1 lame2 thick];

% Tolerance
tol = 1e-9;
max_iters = 50;

% Total / Updated Lagrangian
framework = 'UL';