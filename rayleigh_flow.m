function [m_dot,M_3,p_3,T_stag_3,T_3] = rayleigh_flow(M_2,T_2,p_2,R,hi,width)

m_dot = ((p_2)/(R*T_2))*(M_2)*sqrt(1.4*R*T_2)*hi*65*width;
%eq. ratio
phi = 0.95;
cp = 1003; %air
MW_air = 28.85; MW_fuel = 2.00;
%fuel = 'C0H2';
x = 0; y = 2; a = 0.25*(x+y);
n_stoic = 4.76*a;
m_stoic = n_stoic*(MW_air/MW_fuel);

m_dot_f = (m_dot/m_stoic)*phi;
 
%heat release
LHV = 120e6;
n = 0.85; %efficiency

q_heat = n*LHV*(1/(m_dot+m_dot_f));
 
%static/stagnation values at combustor entry and exit
TR = 1/(1+0.2*(M_2)^2);
T_stag_2 = T_2/TR;
T_stag_3 = (q_heat/cp) + T_stag_2;

%exit mach number 
A = 4.8*(1+0.2*M_2^2)*M_2^2; B = (1+1.4*M_2^2)^2;
TR_crit_2 = A/B;

TR_crit_3 = (T_stag_3/T_stag_2)*(TR_crit_2);

%false position for exit mach number (eq..(1))
xl(1) = 1.01;
xu(1) = M_2;

for i = 1:20
    
    fxl(i) = (TR_crit_3/4.8) - ((1+0.2*xl(i)^2)*xl(i)^2)/(1+1.4*xl(i)^2)^2;
    fxu(i) = (TR_crit_3/4.8) - ((1+0.2*xu(i)^2)*xu(i)^2)/(1+1.4*xu(i)^2)^2;
    xr(i) = (xu(i)*fxl(i) - xl(i)*fxu(i)) / (fxl(i)-fxu(i));
    fxr(i) = (TR_crit_3/4.8) - ((1+0.2*xr(i)^2)*xr(i)^2)/(1+1.4*xr(i)^2)^2;
    
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
M_3 = xr(i);
k2 = i;

p_3 = p_2*(1+1.4*M_2^2)/(1+1.4*M_3^2);

T_3 = T_stag_3*(1+0.2*M_3^2)^(-1);




end