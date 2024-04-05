function [h_vec,T_vec,p_vec,rho_vec,a_vec]=isa_prop(h_G_vec)
g_0=9.80665;
R=287.04;
r=6.356766e6;
gamma=1.4;

h_G0_row=[0,11,25,47,53,79,90,105]*1e3;
N_layers=length(h_G0_row)-1;
a_0_row=[-.0065,.003,-.0045,.004];
T_0_row=nan(1,N_layers);
T_0_row(1)=288.16;
p_0_row=nan(1,N_layers);
p_0_row(1)=101330;

% Fill in T_0_row and p_0_row
for n_layer=1:N_layers-1
    [T_0_row(n_layer+1),p_0_row(n_layer+1)]=graient_isothermal_T_p(h_G0_row(n_layer+1),n_layer);
end
[T_0_row.',p_0_row.'] %#ok<NOPRT> 

h_vec=r.*h_G_vec./(r+h_G_vec);
T_vec=nan(size(h_G_vec));
p_vec=T_vec;

for n=1:numel(h_G_vec)
    if h_G_vec(n)<h_G0_row(1) || h_G_vec(n)>h_G0_row(end)
        warning("Invalid value in the input. Values must lie in the region ["+h_G0_row(1)+","+h_G0_row(end)+"]");
    end
    
    for n_layer=1:N_layers
        if h_G_vec(n)<=h_G0_row(n_layer+1)
            [T_vec(n),p_vec(n)]=graient_isothermal_T_p(h_G_vec(n),n_layer);
            break
        end
    end
end
rho_vec=p_vec./R./T_vec;
a_vec=sqrt(gamma.*R.*T_vec);

function [T,p]=graient_isothermal_T_p(h_G,n_layer)
    if mod(n_layer,2)~=0  %n_layer is odd ==> gradient layer
        T=T_0_row(n_layer)+a_0_row((n_layer+1)/2).*(h_G-h_G0_row(n_layer));
        p=p_0_row(n_layer).*(T./T_0_row(n_layer)).^(-g_0./a_0_row((n_layer+1)/2)./R);
    else    %n_layer is even ==> isothermal layer
        T=T_0_row(n_layer);
        p=p_0_row(n_layer).*exp(-g_0.*(h_G-h_G0_row(n_layer))./R./T_0_row(n_layer));
    end
end
end