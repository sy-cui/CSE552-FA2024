% Solve Partitioned Finite Element Matrix System
%
% Copyright (C) Arif Masud and Tim Truster
% 7/2009
% UIUC

%Move Constrained DOF to RHS
Fdtilda = zeros(neq,1);

for i = 1:neq
    rhs = 0;
    for j = 1:nieq
       rhs = rhs + Kdf(i,j)*ModelDc(j);
    end
    Fdtilda(i) = Fd(i) - rhs;
end

% Compute fint and map it to the equation numbering
FintStressStrain;

rhs = Fdtilda - Fdint;
r0 = norm(rhs); res = 1;
iters = 0;
while res > tol && iters < 100
    % Solve Kd = F
    du = Kdd\rhs;
    ModelDx = ModelDx + du; 
    
    % Update residual
    FintStressStrain;

    % Fdint
    rhs = Fdtilda - Fdint;
    res = norm(rhs) / r0;

    % for node = 1:numnp
    %     for dir = 1:ndf
    %         gDOF = NDOFT(node, dir);
    %         if gDOF <= neq
    %             u(dir, node) = ModelDx(gDOF,1);
    %         else
    %             u(dir, node) = ModelDc(gDOF - neq);
    %         end
    %     end
    % end
    % FormFE;

    iters = iters + 1;
end

disp(['Iters: ' num2str(iters)])
