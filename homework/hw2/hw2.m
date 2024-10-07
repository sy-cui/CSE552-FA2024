clear; close all;

%% Define thermal conductivity and its derivative
function k = kappa(u, ~, k0, eps)
    k = k0 * exp(eps * u);
end

function dk = dkappa(u, ~, k0, eps)
    dk = (eps * k0) * exp(eps * u);
end

%% Element-level operators

% Element-level internal force vector
function Fint = elem_fint(x, u, w, C, D, Je, kfunc)
    y = reshape(x, [], 1); d = reshape(u, [], 1);
    kval = kfunc(C*d, C*y); % isoparametric
    Fint = (1/Je) * D' * ((kval.*w) .* (D*d));
end

% Element-level consistent tangent
function K = elem_ctan(x, u, w, C, D, Je, kfunc, dkfunc)
    y = reshape(x, [], 1); d = reshape(u, [], 1);
    ue = C*d; xe = C*y;
    kval = kfunc(ue, xe); dkval = dkfunc(ue, xe);
    K = D' * diag(w.*kval) * D; % symmetric part
    K = K + D' * diag(w.*dkval.*(D*d)) * C; % non-symmetric part
    K = (1/Je) * K; % scale with inverse Jacobian
end

% Element-level external force vector (unused)
function Fext = elem_fext(f, w, C, Je)
    g = reshape(f, [], 1);
    Fext = Je * C' * (w.*(C * g));
end

%% Global assemblies

% Assembled internal force vector
function Fint = glob_fint(ne, x, u, w, C, D, Je, kfunc)
    y = reshape(x, [], 1); d = reshape(u, [], 1);
    s = fix(length(y) / ne);
    if s * ne + 1 ~= length(d) || s * ne + 1 ~= length(y)
        error("Size mismatch");
    end

    Fint = zeros(length(y), 1);
    
    for e = 1:ne % Iterate through elements
        i1 = (e - 1) * s + 1;
        i2 = e * s + 1;
        Fint(i1:i2) = Fint(i1:i2) + elem_fint( ...
            x(i1:i2),u(i1:i2),w,C,D,Je,kfunc);
    end
end

% Assembled consistent tangent
function K = glob_ctan(ne, x, u, w, C, D, Je, kfunc, dkfunc)
    y = reshape(x, [], 1); d = reshape(u, [], 1);
    s = fix(length(y) / ne);
    if s * ne + 1 ~= length(d) || s * ne + 1 ~= length(y)
        error("Size mismatch");
    end

    K = zeros(length(y), length(y));
    
    for e = 1:ne % Iterate through elements
        i1 = (e - 1) * s + 1;
        i2 = e * s + 1;
        K(i1:i2, i1:i2) = K(i1:i2, i1:i2) + elem_ctan( ...
            x(i1:i2),u(i1:i2),w,C,D,Je,kfunc,dkfunc);
    end
end

% Assembled external force vector (unused)
function Fext = glob_fext(ne, f, w, C, Je)
    g = reshape(f, [], 1);
    s = fix(length(g) / ne);
    if s * ne + 1 ~= length(g)
        error("Size mismatch");
    end

    Fext = zeros(length(g), 1);
    
    for e = 1:ne % Iterate through elements
        i1 = (e - 1) * s + 1;
        i2 = e * s + 1;
        Fext(i1:i2) = Fext(i1:i2) + elem_fext(g(i1:i2),w,C,Je);
    end
end

%% Solution procedures
k0 = 1; eps = 0.1;
kfunc = @(u, x) kappa(u, x, k0, eps);
dkfunc = @(u, x) dkappa(u, x, k0, eps);

% To use quadratic basis with two-point quadrature, simply set
% p = [-1 0 1];
% [z. w] = gl(1);
p = [-1 1]; % linear basis
[z, w] = gl(0); % one-point quadrature
C = interp_mat(z, p);
D = C * deriv_mat(p);

ne = 4; L = 1; s = length(p) - 1;
x = linspace(0, L, s*ne+1);
he = L / ne; Je = he / 2;
F = zeros(s*ne + 1, 1);

I = speye(s*ne+1); R = I(1:end-1,:);

nlf = @(u) R*glob_fint(ne, x, R'*u, w, C, D, Je, kfunc);
ctf = @(u) R*glob_ctan(ne, x, R'*u, w, C, D, Je, kfunc, dkfunc)*R';

h = 1:50;
un = zeros(length(h), s*ne+1);
u = R*zeros(s*ne+1, 1);
for i = 1:length(h)
    F(1) = h(i);
    [u,k] = nr(u,R*F,nlf,ctf,10,1e-10);
    un(i,:) = R'*u;
end

%% Plotting
function y = soln(x, h, k0, eps)
    h = reshape(h, [], 1); x = reshape(x, 1, []);
    if eps == 0
        y = h ./ k0 .* (1 - x); 
    else
        y = (1/eps).*log(1 + (eps/k0)*h.*(1-x));
    end
end
xf = linspace(0, 1, 65);
ue = soln(xf, h, k0, eps);

figure(1)
subplot(1,2,1)
hold on
    plot(x, un(25,:), '.-k', LineWidth=2, MarkerSize=20, ...
        DisplayName='$N=25$, computed')
    plot(x, un(50,:), '.-r', LineWidth=2, MarkerSize=20, ...
        DisplayName='$N=50$, computed')
    plot(xf, ue(25,:), '--b', LineWidth=2, DisplayName='$N=25$, exact')
    plot(xf, ue(50,:), '--g', LineWidth=2, DisplayName='$N=50$, exact')
hold off; box on; grid on;
xlim([0, L]);
xlabel('Position $x$', Interpreter='latex')
ylabel('Temperature $u(x)$', Interpreter='latex')
legend(Location="northeast", Interpreter='latex')
set(gca, Fontsize=20, FontName='Times new roman')

subplot(1,2,2)
hold on
    plot(h, un(:,1), '.-k', LineWidth=2, MarkerSize=10, ...
        DisplayName='$x=0$, computed')
    plot(h, un(:,s*ne/2+1), '.-r', LineWidth=2, MarkerSize=10, ...
        DisplayName='$x=0.5$, computed')
    plot(h, ue(:,1), '--b', LineWidth=2, DisplayName='$x=0$, exact')
    plot(h, ue(:,33), '--g', LineWidth=2, DisplayName='$x=0.5$, exact')
hold off; box on; grid on;
xlim([h(1) h(end)]);
xlabel('Load steps $n$', Interpreter='latex')
ylabel('Temperature $u(x)$', Interpreter='latex')
legend(Location="northwest", Interpreter='latex')
set(gca, Fontsize=20, FontName='Times new roman')
papersize = [1080 360];
pos = get(gcf, 'Position');
set(gcf, PaperUnits='points', Position=[10 10 papersize], ...
    PaperSize=papersize);
% print -dpdf hw2_eps_0.pdf -bestfit
