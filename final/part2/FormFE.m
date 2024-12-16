% Assemble Quantities from Model Routine
%
% Copyright (C) Arif Masud and Tim Truster
%
% 7/2009
% UIUC

Kdd = zeros(neq,neq);
Kdf = zeros(neq,nieq);
Kfd = zeros(nieq,neq);
Kff = zeros(nieq,nieq);

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
    [xc, ~]        = extract_elem_var(NodeCurr,elem,ndm,nel,ix);
    
    % Extract patch material properties
    ma = ix(elem,nen+1);
    mateprop = MateT(ma,:);

    % Compute and Assemble Patch Stiffness
    EDOFT = LocToGlobDOF(ElemFlag, NDOFT, nel, ndf);
    
    switch iel
        case 1 % Small-Deformation Isotropic Elastostatics Element
            if ndm == 3
                L_Elem1_3d
            else % ndm == 2
                [ElemK,ElemF] = tang_elem_lin_2d(xr,mateprop,nel,ndf,4,PSPS);
            end
        case 2 % Small-strain nonlinear elastostatic element (MIDTERM)
            if ndm == 2
                [ElemK,ElemF] = tang_elem_ss_2d(xc,xr,mateprop,nel,ndf,4);
            else % ndm == 2
                
            end
        case 3 % Stabilized Mixed Pressure-Displacement Element
            if ndm == 2
                L_Elem3_2d
            else % ndm == 3

            end
        case 4 % Implicit Error Element
            if ndm == 2
                L_Elem4_2d
            else % ndm == 3

            end
        case 5 % Stabilized Mixed Pressure-Displacement Element, Error
            if ndm == 2
                L_Elem5_2d
            else % ndm == 3

            end

        case 6 % Neo-Hookean Hyperelastic (FINAL PART I)
            if ndm == 2
                [ElemK,ElemF]=tang_elem_nh_2d(xc,xr,mateprop,nel,ndf);
            else % ndm == 3

            end
    end

    % Assemble Element contribution to Model Quantity
    AssemStifForc;
    
end