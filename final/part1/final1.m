clear; close all; 
header;

%% Input file
input_nh_stretch;
% input_nh_shear;

TRACK_ELEM_IDX = 1; % Which element to track for stress and strain
TRACK_QUAD_IDX = 4; % Which quadrature point to track for stress and strain

%% Run load steps 
steps = 100;
P = linspace(0, 1e6, steps);

sigma_11 = zeros(steps,1);
sigma_22 = zeros(steps,1);
sigma_12 = zeros(steps,1);

E_11 = zeros(steps,1);
E_22 = zeros(steps,1);
E_12 = zeros(steps,1);

d1 = zeros(steps,1);
d2 = zeros(steps,1);

for ii = 1:steps
    NodeLoad(:,3) = P(ii);
    FEA_Program;
    stress_strain;

    if ii == 10
        res10 = res;
    elseif ii == 50
        res50 = res;
    end

    sigma_11(ii) = sigma(1,1);
    sigma_22(ii) = sigma(2,2);
    sigma_12(ii) = sigma(1,2);

    E_11(ii) = GLS(1,1);
    E_22(ii) = GLS(2,2);
    E_12(ii) = GLS(1,2);
    
    u = NodeCurr - NodeTable;
    d1(ii) = u(3,1);
    d2(ii) = u(3,2);
end

%% Plotting
tsteps = linspace(0,1,steps);

lw=2;
mks = 20;
fs = 26;

figure(1)

tiledlayout(1,3)
nexttile

hold on
plot(tsteps, sigma_11.*1e-6, '.-k', 'LineWidth', lw, 'MarkerSize', mks, ...
    'DisplayName', '$\sigma_{11}$');
plot(tsteps, sigma_22.*1e-6, '.-r', 'LineWidth', lw, 'MarkerSize', mks, ...
    'DisplayName', '$\sigma_{22}$');
plot(tsteps, sigma_12.*1e-6, '.-b', 'LineWidth', lw, 'MarkerSize', mks, ...
    'DisplayName', '$\sigma_{12}$');
hold off; box on; grid on;

xlabel('Time / Load steps', 'Interpreter', 'latex');
ylabel('Cauchy stress, [MPa]', 'Interpreter','latex');
legend('Interpreter','latex','Location','northwest')
set(gca, 'fontname', 'Times New Roman', 'Fontsize', fs)

nexttile
hold on
plot(E_11, sigma_11.*1e-6, '.-k', 'LineWidth', lw, 'MarkerSize', mks, ...
    'DisplayName', '$\sigma_{11}$ vs $E_{11}$');
plot(E_22, sigma_22.*1e-6, '.-r', 'LineWidth', lw, 'MarkerSize', mks, ...
    'DisplayName', '$\sigma_{22}$ vs $E_{22}$');
plot(E_12, sigma_12.*1e-6, '.-b', 'LineWidth', lw, 'MarkerSize', mks, ...
    'DisplayName', '$\sigma_{12}$ vs $E_{12}$');
hold off; box on; grid on;

xlabel('Green-Lagrange strain', 'Interpreter', 'latex');
ylabel('Cauchy stress, [MPa]', 'Interpreter','latex');
legend('Interpreter','latex','Location','southeast')
set(gca, 'fontname', 'Times New Roman', 'Fontsize', fs)

nexttile
hold on
plot(tsteps, d1, '.-k', 'LineWidth', lw, 'MarkerSize', mks, ...
    'DisplayName', '$d_1$');
plot(tsteps, d2, '.-r', 'LineWidth', lw, 'MarkerSize', mks, ...
    'DisplayName', '$d_2$');
hold off; box on; grid on;

xlabel('Time / Load steps', 'Interpreter', 'latex');
ylabel('Node 3 Displacement [m]', 'Interpreter','latex');
legend('Interpreter','latex','Location','northwest');
set(gca, 'fontname', 'Times New Roman', 'Fontsize', fs);

papersize = [1640 360];
set(gcf, PaperUnits='points', Position=[10 10 papersize], ...
    PaperSize=papersize);
print -dpdf final1_stretch.pdf -bestfit
