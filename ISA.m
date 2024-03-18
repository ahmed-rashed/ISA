clc
clearvars
close all

g_0=9.80665;
R=287.04;
r=6.356766e6;
gamma=1.4;

N_layer=50;
h_G0_row=[0,11,25,47,53,79,90,105]*1e3;
T_0_row=[288.16,216.66,216.66,282.66,282.66,165.66,165.66];
p_0_row=[101330,22632,2488.6,120.44,58.321,1.0094,.10444];
a_0_row=[-.0065,.003,-.0045,.004];

h_G_row=nan(1,7*N_layer);
T_row=nan(1,7*N_layer);
p_row=nan(1,7*N_layer);
for n=1:7
    ind_layer_vec=N_layer*(n-1)+(1:N_layer);
    h_G_row(ind_layer_vec)=linspace(h_G0_row(n),h_G0_row(n+1),N_layer);

    if mod(n,2)~=0
        ind_a_0=(n+1)/2;
        T_row(ind_layer_vec)=T_0_row(n)+a_0_row(ind_a_0).*(h_G_row(ind_layer_vec)-h_G0_row(n));
        p_row(ind_layer_vec)=p_0_row(n).*(T_row(ind_layer_vec)./T_0_row(n)).^(-g_0./a_0_row(ind_a_0)./R);
    else
        T_row(ind_layer_vec)=repelem(T_0_row(n),1,N_layer);
        p_row(ind_layer_vec)=p_0_row(n).*exp(-g_0.*(h_G_row(ind_layer_vec)-h_G0_row(n))./R./T_0_row(n));
    end
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
plot(p_row/1e5,h_G_row./1e3)
xlabel('p (bar)')
ylabel('h_G (km)')

nexttile
plot(rho_row,h_G_row./1e3)
xlabel('rho (kg/m^3)')
ylabel('h_G (km)')