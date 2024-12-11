function ELDOFT = LocToGlobDOF(NFlags, NDOFT, maxnodes, ndir)
%
% Copyright (C) Arif Masud and Tim Truster
%
% maxnodes: number of nodes in the current element
% NFlagsL: (maxnodes x 1), global nodal index for each elemental node
% NDOFT: (total nodes x 2maxdof) Mapping between node and equation indices
% ndir: nodal dof / directions
% ELDOFT: Compact equation to node numbering 

ELDOFT = zeros(1,ndir*maxnodes);
loc=0;

for node = 1:maxnodes
    
    for direction = 1:ndir
        loc = loc+1;
        ELDOFT(loc) = NDOFT(NFlags(node), direction);
    end
    
end