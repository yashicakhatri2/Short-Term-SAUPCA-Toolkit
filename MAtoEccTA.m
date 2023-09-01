function [OUT1, OUT2] = MAtoEccTA(M,e)
M = wrapToPi(M);
flag = 1;
n = 0;
if (M > -pi && M < 0) || M > pi
    E(1) = M - e;
else
    E(1) = M + e;
end

while (flag == 1) || (abs(E(n+1) - E(n)) > 1E-10 && n < 100)
    n = n + 1;
    E(n+1) = E(n) + (-E(n) + e * sin(E(n)) + M) / (1 - e * cos(E(n)));
    flag = 0;  
end
E_final = wrapToPi(E(end));
TA = 2 * atan(sqrt((1+e)/(1-e)) * tan(E_final/2));

OUT1 = E_final;
OUT2 = TA;
end