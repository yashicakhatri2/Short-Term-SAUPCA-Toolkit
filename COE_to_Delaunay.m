function OUT = COE_to_Delaunay(IN,mu,dir)
% Convert Classical Orbital Elements to a Delaunay set and vice-versa.
% =======================================================================
% INPUTS = 
% IN: Input set - COE or Delaunay set
% mu: Gravitational parameter
% dir: 1 - COE to Delaunay set, 0 - Delaunay set to COE
% 
% VARIABLES = 
% a: Semi-major axis 
% e: Eccentricity
% i: Inclination
% O: Right Ascension of the Ascending Node
% w: Arguemnt of periapsis
% M: Mean Anomaly
% l: Delaunay set element 1
% g: Delaunay set element 2
% h: Delaunay set element 3
% L: Delaunay set element 4
% G: Delaunay set element 5
% H: Delaunay set element 6
%
% OUTPUT =
% OUT: Output set - COE or Delaunay set
% =======================================================================

if dir == 1
    a = IN(1); e = IN(2); i = IN(3); O = IN(4); w = IN(5); M = IN(6);
    l = wrapToPi(M);
    g = wrapToPi(w);
    h = wrapToPi(O);
    L = sqrt(mu*a);
    G = L * sqrt(1 - e^2);
    H = G * cos(i);
    OUT = [l;g;h;L;G;H];
else
    l = IN(1); g = IN(2); h = IN(3); L = IN(4); G = IN(5); H = IN(6);
    a = L^2/mu;
    e = sqrt(1-(G/L)^2);
    i = wrapToPi(acos(H/G));
    O = wrapToPi(h);
    w = wrapToPi(g);
    M = wrapToPi(l);
    OUT = [a;e;i;O;w;M];
end

end