function M  = mass_elem_2d(X,rho,t,nel,ndf,lumped)
    % Compute the 2D element mass matrix based on configuration 'X'. 
    % 
    % Inputs:
    %   X:          (2 x nel) nodal position vector.
    % 
    %   rho:        Material density in the corresponding configuration.
    %
    %   t:          Thickness of the 2D structure in the normal direction.
    % 
    %   nel:        Number of nodes in the element.
    % 
    %   ndf:        Maximum dof per node. Typically 2 in 2D. 
    %   
    %   lumped:     Boolean variable denoting whether use the lumped-mass
    %               diagonal matrix
    % 
    % Outputs:
    %   M:          (nel*ndf x nel*ndf) element mass matrix. Should be a
    %               symmetric matrix consisting of (nel x nel) diagonal
    %               blocks
    %
    
    % Initialize Matrix and Vector
    subM = zeros(nel,nel);
    
    % Load Guass Integration Points
    lint = 4;
    
    % Loop over integration points
    for l = 1:lint
    
        if nel == 3 % Triangle elements
            [Wgt,r,s] =  intpntt(l,lint,0);
        else % Quadralateral elements
            [Wgt,r,s] =  intpntq(l,lint,0);
        end
        
        % Shape functions 
        shp = shpl_2d(r,s,nel);

        % Get element Jacobian
        [~,J] = shpg_2d(shp,X,nel);
    
        % Element mass sub-matrix
        subM = subM + Wgt*rho*J*t*(shp(3,:)' * shp(3,:));
    end

    M = kron(subM, eye(ndf));

    if lumped
        rowsum = sum(M,2);
        M = diag(rowsum);
    end
end