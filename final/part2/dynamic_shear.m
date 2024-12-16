clear; close all;
header;

%% Load input
input_shear;

iel = 6; % Element type
ndf = 2; % Nodal DOF
ndm = 2; % Dimension

TRACK_ELEM_IDX = 1;
TRACK_QUAD_IDX = 4;

%% Static part
max_disp = 0.2;
load_steps = 100;
for lstep = 1:load_steps
    NodeBC(end-1:end,3) = max_disp.*(lstep/load_steps);
    [Fd,ModelDc,neq,nieq,NDOFT] = assign_bc_load_data( ...
        numnp, ndf, NodeBC, NodeLoad, BCLIndex);
    nneq = neq + nieq;
    FormFE;
    SolveFE;
end

disp_ini = NodeCurr - NodeTable

%% Dynamic part
% Modify boundary conditions: remove displacements
BCLIndex = [4 0]';

% Reform the boundary condition
[Fd,ModelDc,neq,nieq,NDOFT] = assign_bc_load_data( ...
    numnp, ndf, NodeBC, NodeLoad, BCLIndex);

% Reform operators
FormMass; 
FormFE;
MKdd = Mdd + (1+alpha)*beta*dt^2*Kdd;
MKdf = Mdf + (1+alpha)*beta*dt^2*Kdf;
MKfd = Mfd + (1+alpha)*beta*dt^2*Kfd;
MKff = Mff + (1+alpha)*beta*dt^2*Kff;

% Initialize velocity and acceleration field
NodeVel = NodeCurr*0.0;

Fint;
ModelA = Mdd \ (Fd - Fdint);
NodeAcc = init_acc(ModelA,NDOFT,neq,numnp,ndf);

tsteps = ceil(tend / dt);
ts = [0 (1:tsteps)*dt];
disp_x = zeros(tsteps+1,1); disp_x(1) = NodeCurr(3,1)-NodeTable(3,1);
disp_y = zeros(tsteps+1,1); disp_y(1) = NodeCurr(3,2)-NodeTable(3,2);
sigma_11 = zeros(tsteps+1,1);
sigma_22 = zeros(tsteps+1,1);
sigma_12 = zeros(tsteps+1,1);

ss_nh_2d;
sigma_11(1) = sigma(1,1);
sigma_22(1) = sigma(2,2);
sigma_12(1) = sigma(1,2);

for tt = 1:tsteps
    % Update K and effective mass matrices
    FormFE;
    MKdd = Mdd + (1+alpha)*beta*dt^2*Kdd;
    MKdf = Mdf + (1+alpha)*beta*dt^2*Kdf;
    MKfd = Mfd + (1+alpha)*beta*dt^2*Kfd;
    MKff = Mff + (1+alpha)*beta*dt^2*Kff;

    % Do step
    time_step;

    % Log data
    disp_x(tt+1) = NodeCurr(3,1)-NodeTable(3,1);
    disp_y(tt+1) = NodeCurr(3,2)-NodeTable(3,2);

    ss_nh_2d;
    sigma_11(tt+1) = sigma(1,1);
    sigma_22(tt+1) = sigma(2,2);
    sigma_12(tt+1) = sigma(1,2);

end


%% Plotting
figure(1)
plot(ts,sigma_11,'-k','LineWidth',1,'DisplayName','$\sigma_{11}$')
hold on 
plot(ts,sigma_22,'-r','LineWidth',1,'DisplayName','$\sigma_{22}$')
plot(ts,sigma_12,'-b','LineWidth',1,'DisplayName','$\sigma_{12}$')
hold off
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Displacement [m]', 'Interpreter', 'latex')
legend('Interpreter','latex','Location','northeast')
set(gca, 'Fontsize', 20, 'Fontname', 'Times new roman')
title(sprintf('$\\alpha = %.2f$', alpha), ...
    'Interpreter', 'latex', 'Fontsize', 24)

papersize = [540 360];
set(gcf, PaperUnits='points', Position=[10 10 papersize], ...
    PaperSize=papersize);
fname = sprintf('final2_shear_%i.pdf', round(-alpha*1000));
print(fname, '-dpdf', '-bestfit')