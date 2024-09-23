clear; close all;

%% Consitutive relation and consistent tangent
% This section contains two functions that computes the nonlinear function
% (the constitutive relation) and its derivative (consistent tangent) at a
% given point u. 
% The functions contain all given parameters in the physical problem,
% ensuring that the rest of the algorithm needs only solve the nonlinear
% equation f(u) = b without involving physical parameters. 

function [F] = constitutive(u)
    L_a = 10; L_b = 5; A = 1; E1 = 1e7; E2 = 1e5; eps = 2e-3;
    eps_a = u / L_a; eps_b = u / L_b;
    ea = eps_a < eps; pa = 1 - ea;
    eb = eps_b < eps; pb = 1 - eb;
    sigma_a = (E1 * eps_a).*ea + (E1 * eps + E2 * (eps_a - eps)).*pa;
    sigma_b = (E1 * eps_b).*eb + (E1 * eps + E2 * (eps_b - eps)).*pb;
    F = A .* (sigma_a + sigma_b);
end

function [K] = tangent(u)
    L_a = 10; L_b = 5; A = 1; E1 = 1e7; E2 = 1e5; eps = 2e-3;
    eps_a = u / L_a; eps_b = u / L_b;
    if eps_a <= eps; K_a = E1; else; K_a = E2; end
    if eps_b <= eps; K_b = E1; else; K_b = E2; end
    K = A * (K_a / L_a + K_b / L_b);
end

%% Problem solve
% We solve the nonlinear problems here using both NR and MNR methods. 
F = [2e4 4e4]; tol=1e-12;

% Newton-Raphson iterations
u = 0; nk = F*0.0; nu = F*0.0; nua = [];
t0 = tic();
for i = 1:length(F)
    [u_arr, ~, ~, k] = nr(u,F(i),@constitutive,@tangent,1e6,tol,1);
    u = u_arr(end);
    nk(i) = k; nu(i) = u; nua = [nua; u_arr];
end
t1 = toc(t0);
disp(['Newton iterations take ' num2str(t1) ' seconds'])

%Modified Newton-Raphson iterations
u = 0.0; mk = F*0.0; mu = F*0.0; mua = [];
t0 = tic();
for i = 1:length(F)
    [u_arr, ~, ~, k] = mnr(u,F(i),@constitutive,@tangent,1e6,tol,1);
    u = u_arr(end);
    mk(i) = k; mu(i) = u; mua = [mua; u_arr];
end
t1 = toc(t0);
disp(['Modified Newton iterations take ' num2str(t1) ' seconds'])


%% Plotting
figure(1)

x = linspace(0, 3e-2, 128);
plot(x, constitutive(x), '-k', LineWidth=2, DisplayName='Exact');
hold on
yline(F, '--k', Linewidth=1, HandleVisibility='off')
scatter(nua, constitutive(nua), 120, 'r', 'filled', DisplayName='N-R')
scatter(mua, constitutive(mua), 60, 'b', 'filled', DisplayName='M-N-R')
hold off; box on;
legend(Interpreter='latex', Location='southeast');
xlabel('Displacement $u$ [cm]', Interpreter='latex')
ylabel('Internal force $F^{int}(u)$ [N]', Interpreter='latex')
set(gca, Fontsize=20, FontName='Times new roman')

papersize = [540 360];
pos = get(gcf, 'Position');
set(gcf, PaperUnits='points', Position=[10 10 papersize], ...
    PaperSize=papersize);
print -dpdf hw1_p3 -bestfit
