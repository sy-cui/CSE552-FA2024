% Solve Partitioned Finite Element Matrix System
%
% Copyright (C) Arif Masud and Tim Truster
% 7/2009
% UIUC

% Compute fint at time t
Fint;
Fdint0 = Fdint;

% Predictors
[NodeCurr,NodeVel,NodeAcc] = predict_kinematics( ...
    NodeCurr,NodeVel,NodeAcc,hht,NodeTable,ModelDc,NDOFT,neq);

% Residual (assume external force is constant for now)
Fint;
rhs = Fd - ((1+alpha)*Fdint - alpha*Fdint0);
r0 = norm(rhs); 

% Multi-correctors
res = 1;
iters = 0;
while res(end) > tol && iters < 1
    % Solve Kd = F
    ModelDa = MKdd\rhs;
    [NodeCurr,NodeVel,NodeAcc] = correct_kinematics( ...
        NodeCurr,NodeVel,NodeAcc,hht,NodeTable,ModelDa,ModelDc,NDOFT,neq);
    
    % Update residual
    Fint;
    [ModelA, ~] = extract_model_var(NodeAcc,NDOFT,neq);
    rhs = Fd - ((1+alpha)*Fdint - alpha*Fdint0) - Mdd*ModelA;
    res = [res norm(rhs) / r0]; %#ok<AGROW>

    % Newton-Raphson
    % FormFE;
    % MKdd = Mdd + beta*dt^2*Kdd;
    % MKdf = Mdf + beta*dt^2*Kdf;
    % MKfd = Mfd + beta*dt^2*Kfd;
    % MKff = Mff + beta*dt^2*Kff;

    iters = iters + 1;
end

% disp(['Iters: ' num2str(iters)])
