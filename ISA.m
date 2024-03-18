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


h_G1_row=linspace(h_G0_row(1),h_G0_row(2),N_layer);
h_G2_row=linspace(h_G0_row(2),h_G0_row(3),N_layer);
h_G3_row=linspace(h_G0_row(3),h_G0_row(4),N_layer);
h_G4_row=linspace(h_G0_row(4),h_G0_row(5),N_layer);
h_G5_row=linspace(h_G0_row(5),h_G0_row(6),N_layer);
h_G6_row=linspace(h_G0_row(6),h_G0_row(7),N_layer);
h_G7_row=linspace(h_G0_row(7),h_G0_row(8),N_layer);

h_G_row=[h_G1_row,h_G2_row,h_G3_row,h_G4_row,h_G5_row,h_G6_row,h_G7_row];

h_row=r.*h_G_row./(r+h_G_row);

T1_row=T_0_row(1)+a_0_row(1).*(h_G1_row-h_G0_row(1));
T2_row=repelem(T_0_row(2),1,N_layer);
T3_row=T_0_row(3)+a_0_row(2).*(h_G3_row-h_G0_row(3));
T4_row=repelem(T_0_row(4),1,N_layer);
T5_row=T_0_row(5)+a_0_row(3).*(h_G5_row-h_G0_row(5));
T6_row=repelem(T_0_row(6),1,N_layer);
T7_row=T_0_row(7)+a_0_row(4).*(h_G7_row-h_G0_row(7));

T_row=[T1_row,T2_row,T3_row,T4_row,T5_row,T6_row,T7_row];

p1_row=p_0_row(1).*(T1_row./T_0_row(1)).^(-g_0./a_0_row(1)./R);
p2_row=p_0_row(2).*exp(-g_0.*(h_G2_row-h_G0_row(2))./R./T_0_row(2));
p3_row=p_0_row(3).*(T3_row./T_0_row(3)).^(-g_0./a_0_row(2)./R);
p4_row=p_0_row(4).*exp(-g_0.*(h_G4_row-h_G0_row(4))./R./T_0_row(4));
p5_row=p_0_row(5).*(T5_row./T_0_row(5)).^(-g_0./a_0_row(3)./R);
p6_row=p_0_row(6).*exp(-g_0.*(h_G6_row-h_G0_row(6))./R./T_0_row(6));
p7_row=p_0_row(7).*(T7_row./T_0_row(7)).^(-g_0./a_0_row(4)./R);

p_row=[p1_row,p2_row,p3_row,p4_row,p5_row,p6_row,p7_row];

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