%test various functions with this code

%% CTRLSURF
DTOR = pi/180;
Se = 5*DTOR; alpha = 4.0*DTOR;
RTOD = 180/pi;
xc = -0.70; zc = -0.10;
M_1(1) = 8.0; h(1) = 16000;

for i = 1:100
    
    [p_1 ,rho_1, T_1 ,a] = pvnrt(h);

    [F_xc(i), F_zc(i), Mc(i)] = ctrlsurf(Se,alpha,M_1(i),p_1,RTOD,DTOR,zc,xc);
    
    M_1(i+1) = M_1(i) + 0.05;

end

%% SOMETHING ELSE

