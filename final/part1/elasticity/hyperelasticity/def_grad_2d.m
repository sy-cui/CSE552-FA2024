function [F,J] = def_grad_2d(r,s,x,X)
    % Compute the 2D deformation gradient tensor F, defined as 
    %       F_ij = \partial x_i / \partial X_j,
    % at the reference points (r, s).
    % 
    % Inputs:
    % 
    %   r:          First coordinate in the reference domain.
    % 
    %   s:          Second coordinate in the reference domain.
    %
    %   x:          Current nodal position. Shape should be (ndm x nel)
    %
    %   X:          Reference nodal position. Shape should be (ndm x nel)
    %
    % Outputs:
    % 
    %   F:          (2 x 2) deformation gradient tensor.
    %
    %   J:          Jacobian / determinant of F.
    %

    nel = size(x,2);
    shp = shpl_2d(r,s,nel);
    [Qxy, ~] = shpg_2d(shp,X,nel);
    F = x * Qxy(1:2,:)';
    J = F(1,1)*F(2,2) - F(1,2)*F(2,1);
end