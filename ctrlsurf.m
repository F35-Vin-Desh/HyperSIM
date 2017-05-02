function [F_xc, F_zc, Mc] = ctrlsurf(Se,alpha,M_1,p_1,RTOD,DTOR,zc,xc)
%M_1,p_1,T_1
%solve the turn angle, 
if abs(Se) == abs(alpha) && Se == -1*alpha
    del_P = 0; F_xc = 0; F_zc = 0; M_c = 0; TURN = 0; X = TURN;
end
 
%LET 1 REPRESENT PTOD AND 2 REPRESENT PDOT
set = 0; %this variable sets it
if abs(Se) > abs(alpha)
        
    if alpha >= 0 && Se >= 0
    TURN = abs(Se)+abs(alpha);
    set = 1;

    elseif alpha >= 0 && Se <= 0
    TURN = abs(Se)-abs(alpha);
    set = 2;

    elseif alpha <= 0 && Se >= 0
    TURN = abs(Se) - abs(alpha);
    set = 1;

    else
    TURN = abs(Se)+abs(alpha);
    set = 2;
    
    end
    
end
  
if abs(Se) < abs(alpha)
    if alpha >= 0 && Se >= 0
    TURN = abs(alpha)+abs(Se);
    set = 2;

    elseif alpha >= 0 && Se <= 0
    TURN = abs(alpha)-abs(Se);
    set = 1;

    elseif alpha <= 0 && Se >= 0
    TURN = abs(alpha) - abs(Se);
    set = 2;

    else
    TURN = abs(Se)+abs(alpha);
    set = 1;
    end
end
X = TURN; %in radians already.- ensure LESS than 15 deg.
 
%oblique calculation
b = -(M_1^2+2)/M_1^2 - 1.4*(sin(X))^2;
c = (2*M_1^2+1)/M_1^4 + (1.44 + 1.4/M_1^2)*(sin(X))^2;
d = -((cos(X))^2)/M_1^4;

POLY = [1 b c d]; 
r = roots(POLY);

%find the middle real root
A = zeros(1,3);
for k = 1:length(A)
    %check if real
    if isreal(r(k))==1
        A(k) = asin(sqrt(r(k)))*RTOD;
    end
    %check if makes sense
    if A(k) <= X*RTOD
        A(k) = [];
    elseif (X*RTOD<A(k)) && (A(k)<75)
        A(k) = A(k);
    else
        A(k) = [];
    end
end
THT_CS = max(A)*DTOR;
S = (sin(THT_CS))^2; %THT_CS = oblique shock from ctrl surface, angle
po_ds = (p_1/6)*(7*M_1^2*S-1); %po - pressure due to oblique (O) shock

%expansion wave calculation - given M_1 and p_1
v_2 = X + sqrt(6)*atan(sqrt((M_1^2-1)/6))-atan(sqrt(M_1^2-1)); %v_2 = X + f(M)

%find M_ods from v_2- 
xl(1) = M_1;
xu(1) = 4*M_1;

for k = 1:500
    
    fxl(k) = v_2 - (sqrt(6)*atan(sqrt((xl(k)^2-1)/6))-atan(sqrt(xl(k)^2-1)));
    fxu(k) = v_2 - (sqrt(6)*atan(sqrt((xu(k)^2-1)/6))-atan(sqrt(xu(k)^2-1)));
    xr(k) = (xu(k)*fxl(k) - xl(k)*fxu(k)) / (fxl(k)-fxu(k));
    fxr(k) = v_2 - (sqrt(6)*atan(sqrt((xr(k)^2-1)/6))-atan(sqrt(xr(k)^2-1)));

    if k >=2
        ea(k) = ((xr(k) - xr(k-1))/xr(k))*100;
        if ea(k)<0.005;
            break
        end
    %plot(1:i,ea,'LineWidth',1.5);
    end
    %compare conditions for root finding

    if fxl(k)*fxr(k)<0
        xl(k+1) = xl(k);
        xu(k+1) = xr(k);

    elseif fxl(k)*fxr(k)>0
        xl(k+1) = xr(k);
        xu(k+1) = xu(k);

    else
        break
    end

end
M_pds = xr(k);
pp_ds =  ((1+0.2*M_1^2)/(1+0.2*M_pds^2))^3.5; %pp_ds = pressure downstream due to prandtl-meyer expansion

if set == 1
    pcl = po_ds; 
    pcu = pp_ds;
else
    pcl = pp_ds;
    pcu = po_ds;
end
Ar = 0.25; %control surface area.

F_zc = -(pcl-pcu)*cos(X)*Ar;
F_xc = -(pcl-pcu)*sin(X)*Ar;

Mc = zc*F_xc - xc*F_zc;
end


