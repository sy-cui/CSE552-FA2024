clear; close all;

%% Consitutive relation and consistent tangent
function [F] = constitutive(u)
    F = (0.19*u.^3 - 2*u.^2 + 6*u + 0.1) .* exp(0.02.*u);
end

function [K] = tangent(u)
    K = (0.57*u.^2 - 4*u + 6) * exp(0.02.*u) + 0.02*constitutive(u);
end

%% Nonlinear Newton iterations
function [u, iter, res] = newton(u, F, tol)
    res = F - constitutive(u);
    iter = 0;
    while norm(res) > tol
        K = tangent(u);
        du = K\res;
        u = u + du;
        res = F - constitutive(u);
        iter = iter + 1;
    end
end

function [u, iter, res] = modified_newton(u, F, tol)
    res = F - constitutive(u);
    iter = 0;
    K = tangent(u); % This is the only difference
    while norm(res) > tol
        du = K\res;
        u = u + du;
        res = F - constitutive(u);
        iter = iter + 1;
    end
end

%% Problem solve
F = 0.5:0.5:5.5; tol=1e-12;

% Newton-Raphson iterations
u = 0; nk = F*0.0; nu = F*0.0; nua = []; nra = [];
t0 = tic();
for i = 1:length(F)
    [u_arr, r_arr, ~, k] = nr(u,F(i),@constitutive,@tangent,1e6,tol,1);
    u = u_arr(end);
    nk(i) = k; nu(i) = u; nua = [nua; u_arr];  nra = [nra; r_arr];
end
t1 = toc(t0);
disp(['Newton iterations take ' num2str(t1) ' seconds'])

%Modified Newton-Raphson iterations
u = 0.0; mk = F*0.0; mu = F*0.0; mua = []; mra = [];
t0 = tic();
for i = 1:length(F)
    [u_arr, r_arr, ~, k] = mnr(u,F(i),@constitutive,@tangent,1e6,tol,1);
    u = u_arr(end);
    mk(i) = k; mu(i) = u; mua = [mua; u_arr];  mra = [mra; r_arr];
end
t1 = toc(t0);
disp(['Modified Newton iterations take ' num2str(t1) ' seconds'])

%% Plotting
figure(1)
subplot(1,2,1)
x = linspace(0, 2, 128);
plot(x, constitutive(x), '-k', LineWidth=2, DisplayName='Exact');
hold on
yline(F, '--k', Linewidth=1, HandleVisibility='off')
scatter(nua, constitutive(nua), 120, 'r', 'filled', DisplayName='N-R')
scatter(mua, constitutive(mua), 60, 'b', 'filled', DisplayName='M-N-R')
hold off; box on;
legend(Interpreter='latex', Location='southeast');
xlabel('Displacement $u$', Interpreter='latex')
ylabel('Internal force $F^{int}(u)$', Interpreter='latex')
set(gca, Fontsize=20, FontName='Times new roman')

subplot(1,2,2)
scatter(1:length(F), nk, 40, 'r', 'filled', DisplayName='Newton')
hold on
scatter(1:length(F), mk, 40, 'b', 'filled', DisplayName='Modified Newton')
hold off; box on;
legend(Interpreter='latex', Location='northwest');
xlabel('Load steps', Interpreter='latex')
ylabel('Iterations', Interpreter='latex')
set(gca, Fontsize=20, FontName='Times new roman')

papersize = [1080 360];
pos = get(gcf, 'Position');
set(gcf, PaperUnits='points', Position=[10 10 papersize], ...
    PaperSize=papersize);
print -dpdf hw1_p4 -bestfit

figure(2)
subplot(2,1,1)
semilogy(1:sum(nk)+length(nk), abs(nra), '*-r', LineWidth=1.5,MarkerSize=5)
hold on; yline(tol, '--k', linewidth=1.5); hold off;
xlabel('Newton-Raphson Iterations', Interpreter='latex')
ylabel('Residual [N]', Interpreter='latex')
set(gca, Fontsize=20, FontName='Times new roman')

subplot(2,1,2)
semilogy(1:sum(mk)+length(mk), abs(mra), '*-b', LineWidth=1.5,MarkerSize=5)
hold on; yline(tol, '--k', linewidth=1.5); hold off;
xlabel('Modified Newton-Raphson Iterations', Interpreter='latex')
ylabel('Residual [N]', Interpreter='latex')
set(gca, Fontsize=20, FontName='Times new roman')

papersize = [720 540];
pos = get(gcf, 'Position');
set(gcf, PaperUnits='points', Position=[10 10 papersize], ...
    PaperSize=papersize);
print -dpdf hw1_p4_res -bestfit
