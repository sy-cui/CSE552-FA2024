% Arbitrary data for assistance in defining the mesh
L = 1;
H = 1;
P = linspace(0, 2, 20);

% Mesh Nodal Coordinates
NodeTable = [0 0
             L 0
             0 H
             L H];
numnp = length(NodeTable);

% Mesh Element Connectivities
ix = [1 2 4 3 1];
nen = 4;
numel = 1;

% Mesh Boundary Conditions and Loads
BCLIndex = [6 2]';
NodeBC = [1 1 0
          1 2 0
          2 1 0
          2 2 0
          3 2 0
          4 2 0];
NodeLoad = [3 1 P
            4 1 P];

% Mesh Material Properties
% young = 10e5;
% pois = .25;
% a = pois * young / (1 + pois) / (1 - 2*pois);
% b = young / (1 + pois);
% c = b/10;
a = 40; 
b = 80;
c = 160;

thick = 1;
PSPS = 'n';
MateT = [a b c thick];
% MateT = [young pois thick];

tol = 1e-12;