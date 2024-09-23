function [u_arr, r_arr, f_arr, iter] = nr( ...
    u, ...              % Initial guess 
    F, ...              % RHS load value
    nlf, ...            % Nonlinear function
    ctf, ...            % Consistent tangent function
    max_iter, ...       % Maximum allowable numbe of iterations
    abs_tol, ...        % Absolute tolerance
    rel_tol ...         % Relative tolerance
)
    iter = 0; Fint = nlf(u); res = F - Fint; dim = length(u);
    tol = min(abs_tol, rel_tol * res);
    
    bsize = 4;
    u_arr = zeros(4, dim); u_arr(1, :) = u;
    r_arr = zeros(4, dim); r_arr(1, :) = res;
    f_arr = zeros(4, dim); f_arr(1, :) = Fint;

    while res > tol && iter < max_iter
        K = ctf(u);
        du = K\res;
        u = u + du;
        Fint = nlf(u);
        res = F - Fint;
        iter = iter + 1;

        if iter > bsize
            new_bsize = ceil(1.5 * bsize);
            new_u_arr = zeros(new_bsize, dim);
            new_r_arr = zeros(new_bsize, dim);
            new_f_arr = zeros(new_bsize, dim);
            new_u_arr(1:bsize, :) = u_arr;
            new_r_arr(1:bsize, :) = r_arr;
            new_f_arr(1:bsize, :) = f_arr;
            u_arr = new_u_arr;
            r_arr = new_r_arr;
            f_arr = new_f_arr;
            bsize = new_bsize;
        end

        u_arr(iter, :) = u;
        r_arr(iter, :) = res;
        f_arr(iter, :) = Fint;
    end
    
    u_arr = u_arr(1:iter);
    r_arr = r_arr(1:iter);
    f_arr = f_arr(1:iter);
    
end