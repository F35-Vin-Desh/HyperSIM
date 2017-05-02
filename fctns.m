function [dudt, dwdt, dxdt, dhdt, dqdt, dadt] = fctns(EF_X, EF_Z, EM, Iyy, m, g, alpha, q, u, w)
    
    dudt = EF_X/m - g*sin(alpha) - q*w;
    
    dwdt = EF_Z/m + g*cos(alpha) + q*u;
    
    dxdt = cos(alpha)*u + sin(alpha)*w;
    
    dhdt = -sin(alpha)*u + cos(alpha)*w;
    
    dqdt = EM/Iyy;
    
    dadt = q;
    
end
