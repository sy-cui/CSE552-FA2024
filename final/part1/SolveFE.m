% Solve Partitioned Finite Element Matrix System
%
% Copyright (C) Arif Masud and Tim Truster
% 7/2009
% UIUC

%Move Constrained DOF to RHS
Fdtilde = zeros(neq,1);

for i = 1:neq
    rhs = 0;
    for j = 1:nieq
       rhs = rhs + Kdf(i,j)*ModelDc(j);
    end
    Fdtilde(i) = Fd(i) - rhs;
end

% Compute fint and map it to the equation numbering
fint_2d;

rhs = Fdtilde - Fdint;
r0 = norm(rhs); 
res = 1;
iters = 0;
while res(end) > tol && iters < max_iters
    % Solve Kd = F
    ModelDx = Kdd\rhs;
    du = model_to_global(ModelDx,ModelDc,NDOFT,numnp,ndf,neq);
    NodeCurr = NodeCurr+du;
    
    % Update residual
    fint_2d;
    rhs = Fdtilde - Fdint;
    res = [res norm(rhs) / r0]; %#ok<AGROW>
    
    % Update tangent if using Newton-Raphson
    % FormFE;

    iters = iters + 1;
end

disp(['Iters: ' num2str(iters)])
