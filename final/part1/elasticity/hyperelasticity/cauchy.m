function S = cauchy(F,J,mateprop)
    % Compute the Cauchy stress tensor from the deformation gradient based 
    % on compressible Neo-Hookean hyperelasticity. 
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
    I = eye(dim);

    S = lame1/J*log(J)*I + lame2/J*(F*F'-I);
end