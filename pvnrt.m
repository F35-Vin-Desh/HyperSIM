%pressure, temperature, and density calculator
function [p, rho, T_K, a] = pvnrt(h)
    
    if (11000>h) && (h<25000)
        T = -56.46; %C
        T_K = T + 273.15;
        p = 22.65*exp(1.73-0.000157*h);
        
    elseif h>=25000
        T = -131.21 + 0.00299*h;
        T_K = T + 273.15;
        p = 2.488*((T+273.1)/216.6)^-11.388;
        
    else 
        T = 15.04 - 0.00649*h;
        T_K = T + 273.15;
        p = 101.29*((T+273.1)/288.08)^5.256;
    end

    rho = p / (0.2869*(T+273.1));
    
    a = sqrt(1.4*287*T_K); %speed of sound

end

        