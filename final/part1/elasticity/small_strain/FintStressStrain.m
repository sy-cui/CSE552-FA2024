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
    
    % xl (dim x nel): nodal position
    xl = zeros(ndm, nel);
    ElemFlag = zeros(nel, 1);
    for k = 1:nel
        node = ix(elem,k);
        ElemFlag(k) = node;
        for l = 1:ndm
            xl(l,k) = NodeTable(node,l);
        end
    end

    % ul (dim x nel): nodal displacement
    EDOFT = LocToGlobDOF(ElemFlag, NDOFT, nel, ndf);
    ul = zeros(ndm,nel);
    for i = 1:nel*ndf
        ndof_index = EDOFT(i);
        if(ndof_index<=neq)
            ul(i) = ModelDx(ndof_index);
        else
            ul(i) = ModelDc(ndof_index-neq);
        end
    end
    
    %Extract patch material properties
    ma = ix(elem,nen+1);
    mateprop = MateT(ma,:);
    
    % Set Material Properties
    a = mateprop(1);
    b = mateprop(2);
    c = mateprop(3);
    
    % Load Guass Integration Points
    if nel == 3
        lint = 4;
    else
        lint = 4;
    end
    
    % Initialize Matrix and Vector
    nst = nel*ndf;
    strain = zeros(lint,3);
    stress = zeros(lint,3);
    ul_elem = reshape(ul,ndf*nel,1);
    fintloc = 0;

    % Loop over integration points
    for l = 1:lint
        if nel == 3
            [Wgt,r,s] =  intpntt(l,lint,0);
        else
            [Wgt,r,s] =  intpntq(l,lint,0);
        end

        % Evaluate local basis functions at integration point
        shp = shpl_2d(r,s,nel);

        % Evaluate first derivatives of basis functions at int. point
        [Qxy, Jdet] = shpg_2d(shp,xl,nel);

        % Form B matrix
        if nel == 3
        Bmat = [Qxy(1,1) 0        Qxy(1,2) 0        Qxy(1,3) 0        
                0        Qxy(2,1) 0        Qxy(2,2) 0        Qxy(2,3)
                Qxy(2,1) Qxy(1,1) Qxy(2,2) Qxy(1,2) Qxy(2,3) Qxy(1,3)];
        else
        
        Bmat = [Qxy(1,1) 0        Qxy(1,2) 0        Qxy(1,3) 0        Qxy(1,4) 0 
                0        Qxy(2,1) 0        Qxy(2,2) 0        Qxy(2,3) 0        Qxy(2,4)
                Qxy(2,1) Qxy(1,1) Qxy(2,2) Qxy(1,2) Qxy(2,3) Qxy(1,3) Qxy(2,4) Qxy(1,4)];
        end
        
        strain_temp = Bmat*ul_elem;
        eps11 = strain_temp(1,1); 
        eps22 = strain_temp(2,1);
        eps12 = strain_temp(3,1)/2;
        
        strain(l,1) = eps11;
        strain(l,2) = eps22;
        strain(l,3) = eps12;

        stress(l,1) = a*(eps11+eps22)+b*eps11+c*(eps11^2+eps12^2);
        stress(l,2) = a*(eps11+eps22)+b*eps22+c*(eps12^2+eps22^2);
        stress(l,3) = b*eps12+c*eps12*(eps11+eps22);

        fintloc = fintloc + Wgt*Jdet*Bmat'*[
            stress(l,1); stress(l,2); stress(l,3)];

    end % quad point loop 

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


