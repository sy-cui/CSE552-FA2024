function Fint = fint_elem_2d(r,s,F,S,X)
    % Compute 2D internal force vector of an element at (quadrature) point
    % (r,s). Continuous formula is 
    %       \int_V W_{i,I} F_{i,J} S_{J,I} dV
    % over the reference domain. 
    %
    % Inputs:
    % 
    %   r:          First coordinate in the reference domain.
    % 
    %   s:          Second coordinate in the reference domain.
    %
    %   F:          (2 x 2) deformation gradient.
    %
    %   S:          (2 x 2) PK2 stress tensor.
    %
    %   X:          Reference nodal position. Shape should be (ndm x nel)
    %
    % Outputs:
    % 
    %   Fint:       Column vector of size (ndm*nel) containing the internal
    %               force vector
    %

    nel = size(X,2);
    shp = shpl_2d(r,s,nel);
    [Qxy, ~] = shpg_2d(shp,X,nel);
    Fint = reshape(F * S * Qxy(1:2,:), [], 1);
end

