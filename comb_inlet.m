%combustor inlet properties, Mach, Pressure, and Temperature

function [M_2,p_2,T_2] = comb_inlet(M_ds,p_ds,T_ds, Ad)

xl(1) = M_ds/2;
xu(1) = 0.975*M_ds;
%*******CONTINUITY EQUATION*
for i = 1:1000
    
    fxl(i) = (xl(i))^2 - (M_ds/Ad)^2*((1+0.2*xl(i)^2)/(1+0.2*M_ds^2))^6;
    fxu(i) = (xu(i))^2 - (M_ds/Ad)^2*((1+0.2*xu(i)^2)/(1+0.2*M_ds^2))^6;
    xr(i) = (xu(i)*fxl(i) - xl(i)*fxu(i)) / (fxl(i)-fxu(i));
    fxr(i) = (xr(i))^2 - (M_ds/Ad)^2*((1+0.2*xr(i)^2)/(1+0.2*M_ds^2))^6;
    
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

%output all values
M_2 = xr(i);
k = i;%convergence speed checker
%pressure and temperature values (equations 7 and 8)
b = (1+0.2*(M_ds)^2)/(1+0.2*(M_2)^2); b = b^3.5;
c = (1+0.2*(M_ds)^2/1+0.2*(M_2)^2);
p_2 = p_ds*b;
T_2 = T_ds*c;


end

    
   
    
    
    
    