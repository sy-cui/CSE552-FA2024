clear; close all;

%% Input file
q4a = 1;
if q4a
    input_quad_stretch;
else
    input_quad_shear;
end

%% FEA_Program.m
steps = length(P);

sigma_11 = zeros(steps,1);
sigma_22 = zeros(steps,1);
sigma_12 = zeros(steps,1);

eps_11 = zeros(steps,1);
eps_22 = zeros(steps,1);
eps_12 = zeros(steps,1);

node = 1;

for ii = 1:steps
    NodeLoad(:,3) = P(ii);
    FEA_Program;

    sigma_11(ii) = stress(node,1);
    sigma_22(ii) = stress(node,2);
    sigma_12(ii) = stress(node,3);

    eps_11(ii) = strain(node,1);
    eps_22(ii) = strain(node,2);
    eps_12(ii) = strain(node,3);
end

%% Post processing
function F = nfunc(eps,sig, a,b,c)
    F(1) = a*(eps(1)+eps(2))+b*eps(1)+c*(eps(1)^2) - sig(1);
    F(2) = a*(eps(1)+eps(2))+b*eps(2)+c*(eps(2)^2) - sig(2);
end

figure(1)
if q4a     
    epse = sigma_22;
    for i = 1:steps
        func = @(x)nfunc(x,[0, sigma_22(i)],a,b,c);
        epss = fsolve(func, [0 0]);
        epse(i) = epss(2);
    end
    
    plot(epse, sigma_22, '-k', LineWidth=2, DisplayName='Exact')
    hold on
    scatter(eps_22,sigma_22,60,'filled','or',DisplayName='FEM')
    hold off; box on; grid on
    xlim([0 0.05]);
    xlabel('Strain $\varepsilon$', Interpreter='latex')
    ylabel('Stress $\sigma$', Interpreter='latex')
    legend(Interpreter='latex', Location='northwest')
    set(gca, Fontsize=24, Fontname='Times New Roman')
    papersize = [540 360];
    set(gcf, PaperUnits='points', Position=[10 10 papersize], ...
        PaperSize=papersize);
    print -dpdf midterm_p4_stretch.pdf -bestfit
else
    hold on
    plot(eps_12, sigma_12, '.-k', Linewidth=2, MarkerSize=20, ...
        DisplayName='$\sigma_{12}$')
    plot(eps_12, sigma_11, '.-r', Linewidth=2, MarkerSize=20, ...
        DisplayName='$\sigma_{11}$')
    plot(eps_12, sigma_22, '.-b', Linewidth=2, MarkerSize=20, ...
        DisplayName='$\sigma_{22}$')
    hold off; box on; grid on
    xlim([0 0.05]);
    xlabel('Strain $\varepsilon$', Interpreter='latex')
    ylabel('Stress $\sigma$', Interpreter='latex')
    legend(Interpreter='latex', Location='northwest')
    set(gca, Fontsize=24, Fontname='Times New Roman')
    papersize = [540 360];
    set(gcf, PaperUnits='points', Position=[10 10 papersize], ...
        PaperSize=papersize);
    print -dpdf midterm_p4_shear.pdf -bestfit
end