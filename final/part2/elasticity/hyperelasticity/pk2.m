function S = pk2(F,J,mateprop)
    % Compute second Piola-Kirchhoff stress tensor from the deformation 
    % gradient based on compressible Neo-Hookean hyperelasticity. 
    % 
    % Inputs:
    %
    %   F:          (d x d) deformation gradient tensor. d = 2 or 3
    %
    %   J:          Determinant of F.
    %
    %   mateprop:   Vector containing material properties. Here we require
    %               that the first and second element are the first and
    %               second Lame parameters \lambda and \mu, respectively.
    %

    dim = size(F,1);
    lame1 = mateprop(1);
    lame2 = mateprop(2);

    invC = inv(F'*F); 
    S = lame1 * log(J) * invC + lame2 * (eye(dim) - invC); %#ok<MINV>
end