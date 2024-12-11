Fdint = zeros(neq, 1);

for elem = 1:numel
    % nel: Number of nodes in current element
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
    
    % Extract patch nodal coordinates
    [xr, ElemFlag] = extract_elem_var(NodeTable,elem,ndm,nel,ix);
    [xc, ~]        = extract_elem_var(NodeCurr,elem,ndm,nel,ix);
    
    % Extract patch material properties
    ma = ix(elem,nen+1);
    mateprop = MateT(ma,:);

    % Compute and Assemble Patch Stiffness
    EDOFT = LocToGlobDOF(ElemFlag, NDOFT, nel, ndf);
    
    % Load Guass Integration Points
    lint = 4;
    
    % Initialize local fint vector
    fintloc = 0;

    % Loop over integration points
    for l = 1:lint
        if nel == 3
            [Wgt,r,s] =  intpntt(l,lint,0);
        else
            [Wgt,r,s] =  intpntq(l,lint,0);
        end

        if strcmpi(framework, 'tl')
            % Evaluate first derivatives of basis functions w.r.t reference
            [Qxy, Jfe] = shpg_2d(shpl_2d(r,s,nel),xr,nel);
            
            % Evaluate deformation gradient
            [F, J] = def_grad_2d(r,s,xc,xr);
    
            % Evaluate PK2 stress tensor
            S = pk2(F,J,mateprop);
    
            % Compute integrand at integration point
            f = fint_elem_2d(r,s,F,S,xr);

        elseif strcmpi(framework, 'ul')
            % Evaluate first derivatives of basis functions w.r.t current
            [Qxy, Jfe] = shpg_2d(shpl_2d(r,s,nel),xc,nel);

            % Evaluate deformation gradient
            [F, J] = def_grad_2d(r,s,xc,xr);
            
            % Compute Cauchy stress
            sig = cauchy(F,J,mateprop);

            % Compute integrand at integration point
            f = reshape(sig * Qxy(1:2,:),[],1);

        else
            error("'framework' must be one of 'ul' and 'tl'")
        end   

        fintloc = fintloc + Wgt*Jfe*mateprop(3)*f;

    end % quad point loop 
    
    % After integration, assemble fintloc into Fdint
    locind1 = 0;
    for ie1 = 1:nel
        for l1 = 1:ndf
            locind1 = locind1 + 1;
            grow = EDOFT(locind1);
            if grow <= neq
                Fdint(grow) = Fdint(grow) + fintloc(locind1);
            end 
        end
    end
    
end % element loop