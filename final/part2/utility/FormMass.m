% Assemble Quantities from Model Routine
%
% Copyright (C) Arif Masud and Tim Truster
%
% 7/2009
% UIUC

Mdd = zeros(neq,neq);
Mdf = zeros(neq,nieq);
Mfd = zeros(nieq,neq);
Mff = zeros(nieq,nieq);

for elem = 1:numel
    
    % Determine element size parameters
    if nen == 3
        nel = 3;
    elseif nen == 4
        if ix(elem,nen) == 0
            nel = 3;
        else
            nel = 4;
        end
    elseif nen == 6
        nel = 6;
    else
        if ix(elem,nen) == 0
            nel = 6;
        else
            nel = 9;
        end
    end
    nst = nel*ndf;
    
    % Extract patch nodal coordinates
    [xr, ElemFlag] = extract_elem_var(NodeTable,elem,ndm,nel,ix);
    
    % Extract patch material properties
    ma = ix(elem,nen+1);
    mateprop = MateT(ma,:);

    % Compute and Assemble Patch Stiffness
    EDOFT = LocToGlobDOF(ElemFlag, NDOFT, nel, ndf);

    % Compute element mass matrix
    ElemM = mass_elem_2d(xr,mateprop(4),mateprop(3),nel,ndf,lumped);

    % Assemble Element contribution to Model Quantity
    AssemMass;
    
end