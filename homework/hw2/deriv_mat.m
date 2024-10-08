function [D] = deriv_mat(p)
    % Compute the differentiation matrix D for Lagrangian nodes p

    n = length(p);
    xi = reshape(p, n, 1);
    diff_ij = (xi - xi') + eye(n);

    diff_ij_inv = 1 ./ diff_ij;
    D_diag_m1 = sum(diff_ij_inv, 2) - 2;

    s = 4 / (max(p) - min(p));
    alpha_i = prod(diff_ij .* s, 2);   
    alpha_i_inv = (1 ./ alpha_i); 

    D = (alpha_i * alpha_i_inv') ./ diff_ij + diag(D_diag_m1);