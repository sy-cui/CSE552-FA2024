% Newton-Raphson Iteration
function [u, k] = nr( ...
    u, ...              % Initial guess 
    F, ...              % RHS load value
    nlf, ...            % Nonlinear function
    ctf, ...            % Consistent tangent function
    max_iter, ...       % Maximum allowable numbe of iterations
    rel_tol ...         % Relative tolerance
)
    k = 0;
    Fint = nlf(u); r = F - Fint; 
    res = norm(r); tol = rel_tol * res;

    while res > tol && k < max_iter
        K = ctf(u);
        du = K\r;
        u = u + du;
        Fint = nlf(u);
        r = F - Fint;
        res = norm(r);
        k = k + 1;
    end
end