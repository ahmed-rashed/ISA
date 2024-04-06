clc
clearvars
close all

h_G_vec=linspace(0,105e3,1e4);

tic; %To measure the execution time
[h_vec,T_vec,p_vec,rho_vec,a_vec]=isa_prop(h_G_vec);
toc %To measure the execution time

figure;
plot(h_vec./1e3,h_G_vec./1e3)
xlabel('$h$ (km)','interpreter','latex')
ylabel('$h_{G}$ (km)','interpreter','latex')

figure;
tiledlayout(1,2)
nexttile
plot(T_vec-273,h_G_vec./1e3)
xlabel('$T$ (\r{ }C)','interpreter','latex')
ylabel('$h_{G}$ (km)','interpreter','latex')

nexttile
plot(a_vec./1e3/(1/60/60),h_G_vec./1e3)
xlabel('$a$ (km/h)','interpreter','latex')
ylabel('$h_{G}$ (km)','interpreter','latex')

figure;
tiledlayout(1,2)
nexttile
plot(p_vec./1e5,h_G_vec./1e3)
xlabel('$p$ (bar)','interpreter','latex')
ylabel('$h_{G}$ (km)','interpreter','latex')

nexttile
plot(rho_vec,h_G_vec./1e3)
xlabel('$\rho$ (kg/m\textsuperscript{3})','interpreter','latex')
ylabel('$h_{G}$ (km)','interpreter','latex')