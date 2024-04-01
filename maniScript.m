clc
clearvars
close all

h_G_vec=linspace(0,105e3,1e4);

tic; %To measure the execution time
[h_vec,T_vec,p_vec,rho_vec,a_vec]=isa_prop(h_G_vec);
toc %To measure the execution time

figure;
plot(h_vec,h_G_vec)
xlabel('$h$ (km)','interpreter','latex')
ylabel('$h_{G}$ (km)','interpreter','latex')

figure;
tiledlayout(1,2)
nexttile
plot(T_vec,h_G_vec)
xlabel('$T$ (\r{ }C)','interpreter','latex')
ylabel('$h_{G}$ (km)','interpreter','latex')

nexttile
plot(a_vec,h_G_vec)
xlabel('$a$ (km/h)','interpreter','latex')
ylabel('$h_{G}$ (km)','interpreter','latex')

figure;
tiledlayout(1,2)
nexttile
plot(p_vec,h_G_vec)
xlabel('$p$ (bar)','interpreter','latex')
ylabel('$h_{G}$ (km)','interpreter','latex')

nexttile
plot(rho_vec,h_G_vec)
xlabel('$\rho$ (kg/m\textsuperscript{3})','interpreter','latex')
ylabel('$h_{G}$ (km)','interpreter','latex')