% Compute internal force 
Fdint = zeros(neq, 1);

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
    [xc, ~]        = extract_elem_var(NodeCurr, elem,ndm,nel,ix);
    
    % Extract patch material properties
    ma = ix(elem,nen+1);
    mateprop = MateT(ma,:);

    % Compute and Assemble Patch Stiffness
    EDOFT = LocToGlobDOF(ElemFlag, NDOFT, nel, ndf);
    
    switch iel
        case 1 % Small-Deformation Isotropic Elastostatics Element
            if ndm == 2
                fint_elem = fint_elem_lin_2d(xc,xr,mateprop,nel,4,PSPS);
            else % ndm == 3
                
            end
        case 2 % Small-strain nonlinear elastostatic element (MIDTERM)
            if ndm == 2
                fint_elem = fint_elem_ss_2d(xc,xr,mateprop,nel,4);
            else % ndm == 2
                
            end
        case 3 % Stabilized Mixed Pressure-Displacement Element
            if ndm == 2
               
            else % ndm == 3
                
            end
        case 4 % Implicit Error Element
            if ndm == 2
                
            else % ndm == 3
                
            end
        case 5 % Stabilized Mixed Pressure-Displacement Element, Error
            if ndm == 2

            else % ndm == 3

            end
        case 6 % Neo-Hookean Hyperelastic (FINAL PART I)
            if ndm == 2
                fint_elem = fint_elem_nh_2d(xc,xr,mateprop,nel,4,framework);
            else % ndm == 3

            end
    end % switch case

    % Assemble fintloc into Fdint
    locind1 = 0;
    for ie1 = 1:nel
        for l1 = 1:ndf
            locind1 = locind1 + 1;
            grow = EDOFT(locind1);
            if grow <= neq
                Fdint(grow) = Fdint(grow) + fint_elem(locind1);
            end 
        end
    end
end