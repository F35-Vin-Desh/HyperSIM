function [M_e,p_e,T_e,THR] = exp_and_thrust(M_3,T_3,p_3,Ae,Ad,An,R,b,w,m_dot,V_1,p_1,p_ds)
%Note: b is a constant defined in 'main'
xl(1) = 1.01*M_3;
xu(1) = An*M_3;

for i = 1:1000
    
    fxl(i) = (xl(i))^2 - (M_3/An)^2*((1+0.2*xl(i)^2)/(1+0.2*M_3^2))^6;
    fxu(i) = (xu(i))^2 - (M_3/An)^2*((1+0.2*xu(i)^2)/(1+0.2*M_3^2))^6;
    xr(i) = (xu(i)*fxl(i) - xl(i)*fxu(i)) / (fxl(i)-fxu(i));
    fxr(i) = (xr(i))^2 - (M_3/An)^2*((1+0.2*xr(i)^2)/(1+0.2*M_3^2))^6;
    
    if i >=2
        ea(i) = ((xr(i) - xr(i-1))/xr(i))*100;
        if ea(i)<0.005;
            break
        end
    %plot(1:i,ea,'LineWidth',1.5);
    end
    %compare conditions for root finding
   
    if fxl(i)*fxr(i)<0
        xl(i+1) = xl(i);
        xu(i+1) = xr(i);
        
    elseif fxl(i)*fxr(i)>0
        xl(i+1) = xr(i);
        xu(i+1) = xu(i);
        
    else
        break
    end

end

M_e = xr(i); k = i;

v = (1+0.2*(M_3)^2)/(1+0.2*(M_e)^2); v = v^3.5;
h = (1+0.2*(M_3)^2)/(1+0.2*(M_e)^2);
p_e = p_3*v;
T_e = T_3*h;

Ve = M_e*sqrt(T_e*1.4*R);
THR = m_dot*(Ve-V_1) + (p_e-p_1)*(Ae/b) - (p_ds-p_1)*(Ae*Ad/b*An);

end
