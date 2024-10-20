clear; close all;

%% Parameters 
tol = 1e-12;    % Relative tolerance
Fext = 3;       % External load
lamb = 0.0;     % Initial load parameter
b = 0.7;        % Arc-length control
da = 0.2;       % Step-size
kmax = 100;     % Maximum MNR iterations
load_steps=35;  % Total number of solution steps

%% Global initalization
u = 0;          % Dispalcement
l = lamb;       % Load parameter
R = 0;          % Main residual 
Du=0; du=0;     % Incremental displacements
Dl=0; dl=0;     % Incremental load parameters
K0 = dfint(u);  % Tangent
q0 = K0 \ Fext; % K\Fext
res = 0;        % Main residual
resa = 0;       % Arc-length residual

c = (1-b) / (q0'*Fext);
K = K0;
Dl_init = da;

%% Modified-Newton-Raphson iterations
data = cell(load_steps,1); % Store data as cell list of structs

for i=1:load_steps
    disp(['Load step ', num2str(i)]);
    cap = 8;
    u_arr = zeros(1, cap);  % displacement 
    l_arr = zeros(1, cap);  % load parameter
    r_arr = zeros(1, cap);  % main residual
    a_arr = zeros(1, cap);  % arc-length residual
    
    % First stored value is the last converged value
    u_arr(1) = u; l_arr(1) = l;
    r_arr(1) = res; a_arr(1) = resa;
    
    % Initialize Dl, Du based on sign of det(K)
    old_K = K; K = dfint(u); q = K \ Fext;
    if det(K) * det(old_K) <= 0
        Dl_init = -Dl_init;    
    end
    Dl = Dl_init; l = l + Dl;
    Du = da*q0;  u = u + Du;

    % Initialize residuals
    R = l*Fext - fint(u); 
    w = K\R;
    r0 = sqrt(abs(R'*w));
    res = 1;
    resa = 0;
    k = 1;
    
    % Second stored value is the initialized u and lambda
    u_arr(2) = u; l_arr(2) = l;
    r_arr(2) = res; a_arr(2) = resa;

    while (resa > 1e-5 || res > tol) && k < kmax
        % Compute incremental u and lambda 
        [dfDu, dfDl] = darc(Du,Dl,K0,b,c);
        dl = (da-arc(Du,Dl,K0,b,c)-dfDu'*w) / (dfDl+dfDu'*q);
        du = w + q*dl;

        % tmp = [K, -Fext; dfDu', dfDl] \ [R; da - arc(Du,Dl,K0,b,c)];
        % du = tmp(1); dl = tmp(2);
        
        % Stability check
        if isnan(du)
            error(['Unstable at load step ', num2str(i)])
        end
        
        % Iteration update
        u = u + du; Du = Du + du;
        l = l + dl; Dl = Dl + dl;
        
        % Update residuals using new values
        R = l*Fext - fint(u);
        w = K\R;
        res = sqrt(abs(R'*w)) / r0;
        resa = abs(arc(Du,Dl,K0,b,c) / da - 1);

        k = k + 1;
        
        % Save data
        if k+1 > cap
            new_cap = ceil(1.7 * cap);

            new_u_arr = zeros(1, new_cap);
            new_l_arr = zeros(1, new_cap);
            new_r_arr = zeros(1, new_cap);
            new_a_arr = zeros(1, new_cap);
            
            new_u_arr(1, 1:cap) = u_arr;
            new_l_arr(1, 1:cap) = l_arr;
            new_r_arr(1, 1:cap) = r_arr;
            new_a_arr(1, 1:cap) = a_arr;
            
            u_arr = new_u_arr;
            l_arr = new_l_arr;
            r_arr = new_r_arr;
            a_arr = new_a_arr;

            cap = new_cap;
        end
        u_arr(k+1) = u;
        l_arr(k+1) = l;
        r_arr(k+1) = res;
        a_arr(k+1) = resa;
    end

    % Remove trailing zeros
    u_arr = u_arr(1:k+1);
    l_arr = l_arr(1:k+1);
    r_arr = r_arr(1:k+1);
    a_arr = a_arr(1:k+1);
    
    % Store data as structure array
    data{i} = struct( ...
        'Iter', k, ...
        'u', u_arr, ...
        'l', l_arr, ...
        'r', r_arr, ...
        'a', a_arr);
end

%% Plotting
t = linspace(-0.5*pi, 0.5*pi, 256);
uc = zeros(length(data), 1);
fc = uc;

figure(1); hold on
uf = linspace(0, 6, 1024);
ff = fint(uf);
plot(uf, ff, '--k', LineWidth=2) % Exact nonlinear function

for i = 1:load_steps
    u = data{i}.u; l = data{i}.l;
    
    % Converged solution
    uc(i) = u(end);
    fc(i) = l(end) * Fext;
    
    % Arc-length function ellipse
    x = u(1) + da / sqrt(c*K0) * cos(t);
    y = (l(1) + da / sqrt(b) * sin(t)) * Fext;
    
    % Solution path
    ud = [u(1) repelem(u(2:end),2)]';
    fd = [reshape([fint(u(1:end-1))', Fext.*l(2:end)']',[],1); 
        fint(u(end))];
    
    plot(x, y, '-g', LineWidth=1)
    plot(ud, fd, '.-b', LineWidth=1.5, MarkerSize=10)
end

scatter(uc, fc, 40, 'red', 'filled')
hold off; box on; grid on
xlim([0, 6]); ylim([-1, 6.5]);
xlabel('$u$', Interpreter='latex')
ylabel('$F$', Interpreter='latex')
set(gca, FontSize=20, FontName='Times New Roman')

papersize = [540 360];
set(gcf, PaperUnits='points', Position=[10 10 papersize], ...
    PaperSize=papersize);
% print -dpdf hw3_b0pt7_soln.pdf -bestfit

%% Plot residual
r = []; a = [];
for i = 1:load_steps
    r = [r, data{i}.r(3:end)];
    a = [a, data{i}.a(3:end)];
end
nr = 1:length(r);
up = round(length(r), -1);

figure(2)
subplot(2,1,1)
semilogy(nr, r, 'k.-', LineWidth=1, MarkerSize=10)
xlim([0, up]);
xlabel('Iterations $k$', Interpreter='latex')
ylabel('$\bar{R}_n^{(k)}$', Interpreter='latex')
set(gca, FontSize=20, FontName='Times New Roman')

subplot(2,1,2)
semilogy(nr, a, 'r.-', LineWidth=1, MarkerSize=10)
xlim([0, up]);
xlabel('Iterations $k$', Interpreter='latex')
ylabel('$\bar{r}_n^{(k)}$', Interpreter='latex')
set(gca, FontSize=20, FontName='Times New Roman')

papersize = [1080 360];
set(gcf, PaperUnits='points', Position=[10 10 papersize], ...
    PaperSize=papersize);
% print -dpdf hw3_b0pt7_res.pdf -bestfit