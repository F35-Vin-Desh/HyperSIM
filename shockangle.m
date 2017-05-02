function [THT_S, DEL] = shockangle(M,T1L, alpha, RTOD)

    %solve the polynomial
    DEL = alpha+T1L;
    b = -(M^2+2)/M^2 - 1.4*(sin(DEL))^2;
    c = (2*M^2+1)/M^4 + (1.44 + 1.4/M^2)*(sin(DEL))^2;
    d = -((cos(DEL))^2)/M^4;
    
    POLY = [1 b c d]; 
    r = roots(POLY);
    
    %find the middle real root
    A = zeros(1,3);
    for i = 1:length(A)
        %check if real
        if isreal(r(i))==1
            A(i) = asin(sqrt(r(i)))*RTOD;
        end
        %check if makes sense
        if A(i) <= DEL*RTOD
            A(i) = [];
        elseif (DEL*RTOD<A(i)) && (A(i)<55)
            A(i) = A(i);
        else
            A(i) = [];
        end
    end
    THT_S = max(A);
end




