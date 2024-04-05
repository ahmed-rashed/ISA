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
for n=1:N_layers-1
    if mod(n,2)~=0  %n is odd ==> gradient layer
        [T_0_row(n+1),p_0_row(n+1)]=graient_T_p(h_G0_row(n+1),h_G0_row(n),T_0_row(n),p_0_row(n),a_0_row((n+1)/2));
    else    %n is even ==> isothermal layer
        [T_0_row(n+1),p_0_row(n+1)]=isothermal_T_p(h_G0_row(n+1),h_G0_row(n),T_0_row(n),p_0_row(n));
    end
end
% [T_0_row.',p_0_row.']

h_vec=r.*h_G_vec./(r+h_G_vec);
T_vec=nan(size(h_G_vec));
p_vec=T_vec;

for n=1:numel(h_G_vec)
    if h_G_vec(n)<h_G0_row(1) || h_G_vec(n)>h_G0_row(end)
        warning("Invalid value in the input. Values must lie in the region ["+h_G0_row(1)+","+h_G0_row(end)+"]");
    end
    
    for m=2:8
        if h_G_vec(n)<=h_G0_row(m)
            if mod(m,2)==0  %m is even ==> gradient layer
                [T_vec(n),p_vec(n)]=graient_T_p(h_G_vec(n),h_G0_row(m-1),T_0_row(m-1),p_0_row(m-1),a_0_row(m/2));
            else    %m is odd ==> isothermal layer
                [T_vec(n),p_vec(n)]=isothermal_T_p(h_G_vec(n),h_G0_row(m-1),T_0_row(m-1),p_0_row(m-1));
            end

            break
        end
    end
end
rho_vec=p_vec./R./T_vec;
a_vec=sqrt(gamma.*R.*T_vec);

    function [T,p]=graient_T_p(h_G,h_G0,T_0,p_0,a_0)
        T=T_0+a_0.*(h_G-h_G0);
        p=p_0.*(T./T_0).^(-g_0./a_0./R);
    end
    
    function [T,p]=isothermal_T_p(h_G,h_G0,T_0,p_0)
        T=T_0;
        p=p_0.*exp(-g_0.*(h_G-h_G0)./R./T_0);
    end
end