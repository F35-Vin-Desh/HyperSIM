%downstream pressure of the oblique shock
%USE ONLY WHEN M,P,T UPSTREAM ARE KNOWN

function [p_ds T_ds M_ds] = oblique_down(p,T_1,M, THT_S, DEL)
   
    %pressure and temperature
    S = (sin(THT_S))^2;
    p_ds = (p/6)*(7*M^2*S-1);
    B = 36*M^2*S; D = (7*M^2*S-1); C = (M^2*S+5);
    T_ds = T_1*(D*C/B);
    
    %mach number downstream
    M_ds = sqrt((1/sin(THT_S-DEL)^2)*(M^2*S+5)/(7*M^2*S-1));
    
end


