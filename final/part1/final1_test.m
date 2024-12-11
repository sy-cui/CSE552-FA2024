clear; close all; 
header;

%% Input file
input_nh_test;

TRACK_ELEM_IDX = 1; % Which element to track for stress and strain
TRACK_QUAD_IDX = 3; % Which quadrature point to track for stress and strain

%% Run load steps 
steps = 100;
P = linspace(0, 1e8, steps);

sigma_11 = zeros(steps,1);
sigma_22 = zeros(steps,1);
sigma_12 = zeros(steps,1);

E_11 = zeros(steps,1);
E_22 = zeros(steps,1);
E_12 = zeros(steps,1);

stretch1 = zeros(steps,1);
stretch2 = zeros(steps,1);

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
    
    D = sqrt(eig(2*GLS+eye(2)));
    stretch1(ii) = D(1);
    stretch2(ii) = D(2);
end

%% Fsolve
function z = nfunc(f12,lame1,lame2,P)
    F1 = f12(1);
    F2 = f12(2);
    J = F1 * F2;
    z(1) = lame1/J*log(J) + lame2/J*(F1^2 - 1);
    z(2) = lame1/J*log(J) + lame2/J*(F2^2 - 1) - 2*P/F1;
end
fsolve_opt = optimset('Display','off');
fsolve_soln = zeros(steps,1);
for ii = 1:steps
    soln = fsolve(@(x) nfunc(x,lame1,lame2,P(ii)), [1; 1], fsolve_opt);
    fsolve_soln(ii) = soln(2);
end

%% Plotting
figure(1)
hold on
plot(fsolve_soln, sigma_22.*1e-6, ...
    'r', 'LineWidth', 2, 'DisplayName', 'Exact');
scatter(stretch2(1:4:end), sigma_22(1:4:end).*1e-6, ...
    60,'k', 'LineWidth', 1.5, 'DisplayName', 'FEM');
hold off; box on; grid on;

xlabel('Principal stretch $\lambda_2$', 'Interpreter', 'latex');
ylabel('Axial true stress $\sigma_{22}$, [MPa]', 'Interpreter','latex');
legend('Interpreter','latex','Location','northwest')
set(gca, 'fontname', 'Times New Roman', 'Fontsize', 20)

set(gca, Fontsize=22, Fontname='Times New Roman')
papersize = [540 360];
set(gcf, PaperUnits='points', Position=[10 10 papersize], ...
    PaperSize=papersize);
print -dpdf final1_test.pdf -bestfit
