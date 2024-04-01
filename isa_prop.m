function [h_vec,T_vec,p_vec,rho_vec,a_vec]=isa_prop(h_G_vec)
g_0=9.80665;
R=287.04;
r=6.356766e6;
gamma=1.4;

h_G0_row=[0,11,25,47,53,79,90,105]*1e3;
T_0_row=[288.16,216.66,216.66,282.66,282.66,165.66,165.66];
p_0_row=[101330,22632,2488.6,120.44,58.321,1.0094,.10444];
a_0_row=[-.0065,.003,-.0045,.004];

h_vec=r.*h_G_vec./(r+h_G_vec);
T_vec=nan(size(h_G_vec));
p_vec=T_vec;

for n=1:numel(h_G_vec)
    if h_G_vec(n)>=h_G0_row(1) && h_G_vec(n)<=h_G0_row(2)
        [T_vec(n),p_vec(n)]=graient_T_p(h_G_vec(n),h_G0_row(1),T_0_row(1),p_0_row(1),a_0_row(1));
    elseif h_G_vec(n)<=h_G0_row(3)
        [T_vec(n),p_vec(n)]=isothermal_T_p(h_G_vec(n),h_G0_row(2),T_0_row(2),p_0_row(2));
    elseif h_G_vec(n)<=h_G0_row(4)
        [T_vec(n),p_vec(n)]=graient_T_p(h_G_vec(n),h_G0_row(3),T_0_row(3),p_0_row(3),a_0_row(2));
    elseif h_G_vec(n)<=h_G0_row(5)
        [T_vec(n),p_vec(n)]=isothermal_T_p(h_G_vec(n),h_G0_row(4),T_0_row(4),p_0_row(4));
    elseif h_G_vec(n)<=h_G0_row(6)
        [T_vec(n),p_vec(n)]=graient_T_p(h_G_vec(n),h_G0_row(5),T_0_row(5),p_0_row(5),a_0_row(3));
    elseif h_G_vec(n)<=h_G0_row(7)
        [T_vec(n),p_vec(n)]=isothermal_T_p(h_G_vec(n),h_G0_row(6),T_0_row(6),p_0_row(6));
    elseif h_G_vec(n)<=h_G0_row(8)
        [T_vec(n),p_vec(n)]=graient_T_p(h_G_vec(n),h_G0_row(7),T_0_row(7),p_0_row(7),a_0_row(4));
    else
        error("Invalid value in the input. Values must lie in the region ["+h_G0_row(1)+","+h_G0_row(end)+"]");
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