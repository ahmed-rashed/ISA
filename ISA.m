clc
clearvars
close all

r=6.356766e6;
gamma=1.4;
R=287.04;
N_layer=50;

h_G0_row=[0,11,25,47,53,79,90,105]*1e3;
N_layers=length(h_G0_row)-1;
T_0_row=nan(1,N_layers);
T_0_row(1)=288.16;
p_0_row=nan(1,N_layers);
p_0_row(1)=101330;

% Fill in T_0_row and p_0_row
for n_layer=1:N_layers-1
    [T_0_row(n_layer+1),p_0_row(n_layer+1)]=graient_isothermal_T_p(h_G0_row,T_0_row,p_0_row,R,h_G0_row(n_layer+1),n_layer);
end
[T_0_row.',p_0_row.'] %#ok<NOPTS> 

h_G_row=nan(1,N_layers*N_layer);
T_row=nan(1,N_layers*N_layer);
p_row=nan(1,N_layers*N_layer);
for n_layer=1:N_layers
    ind_layer_vec=N_layer*(n_layer-1)+(1:N_layer);
    h_G_row(ind_layer_vec)=linspace(h_G0_row(n_layer),h_G0_row(n_layer+1),N_layer);
    [T_row(ind_layer_vec),p_row(ind_layer_vec)]=graient_isothermal_T_p(h_G0_row,T_0_row,p_0_row,R,h_G_row(ind_layer_vec),n_layer);
end
h_row=r.*h_G_row./(r+h_G_row);
rho_row=p_row./R./T_row;
a_row=sqrt(gamma.*R.*T_row);

plot(h_row/1e3,h_G_row./1e3)
xlabel('h (km)')
ylabel('h_G (km)')

figure
tiledlayout(1,2)

nexttile
plot(T_row-273,h_G_row./1e3)
xlabel('T (C)')
ylabel('h_G (km)')

nexttile
plot(a_row./1e3/(1/60/60),h_G_row./1e3)
xlabel('a (km/h)')
ylabel('h_G (km)')

figure
tiledlayout(1,2)

nexttile
plot(p_row./1e5,h_G_row./1e3)
xlabel('p (bar)')
ylabel('h_G (km)')

nexttile
plot(rho_row,h_G_row./1e3)
xlabel('rho (kg/m^3)')
ylabel('h_G (km)')

function [T,p]=graient_isothermal_T_p(h_G0_row,T_0_row,p_0_row,R,h_G_vec,n_layer)
    g_0=9.80665;
    a_0_row=[-.0065,.003,-.0045,.004];

    if mod(n_layer,2)~=0  %n_layer is odd ==> gradient layer
        T=T_0_row(n_layer)+a_0_row((n_layer+1)/2).*(h_G_vec-h_G0_row(n_layer));
        p=p_0_row(n_layer).*(T./T_0_row(n_layer)).^(-g_0./a_0_row((n_layer+1)/2)./R);
    else    %n_layer is even ==> isothermal layer
        T=T_0_row(n_layer);
        p=p_0_row(n_layer).*exp(-g_0.*(h_G_vec-h_G0_row(n_layer))./R./T_0_row(n_layer));
    end
end