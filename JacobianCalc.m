function J = JacobianCalc(COE, mu, t, dir)
% This function calculates the Jacobian between the Equinoctial element set
% and the PQW frame.
%
% =======================================================================
% INPUTS = 
% COE: Current classical orbital element set
% mu: Gravitational parameter of the Earth
% t: Time at which to compute Jacobian
% dir: 1 - Partials of PQW wrt Equinoctial element set, 0 - Partials Equinoctial element set wrt PQW
% 
% VARIABLES = 
% a: Semi-major axis  
% e: Eccentricity
% i: Inclination
% O: Right Ascension of the Ascending Node
% w: Arguemnt of periapsis
% M_t: Mean Anomaly at time, t
% P: Unit vector 1 - [1; 0; 0]
% Q: Unit vector 2 - [0; 1; 0]
% R: Unit vector 3 - [0; 0; 1]
% E: Eccentric anomaly
% r: Distance from Earth
% n: Mean motion
% X: X-term to get PQW state
% Y: Y-term to get PQW state
% Xdot: Xdot-term to get PQW state
% Ydot: Ydot-term to get PQW statedot
% x: State in PQW
% xdot: StateDot in PQW
% L: Transitionary term for partial calcs
% Ldot: Transitionary term for partial calcs
% M: Transitionary term for partial calcs
% Mdot: Transitionary term for partial calcs
% b1: Transitionary term for partial calcs
% b2: Transitionary term for partial calcs
% b3: Transitionary term for partial calcs
% b4: Transitionary term for partial calcs 
% b5: Transitionary term for partial calcs
% b6: Transitionary term for partial calcs
% b1dot: Transitionary term for partial calcs
% b2dot: Transitionary term for partial calcs
% b3dot: Transitionary term for partial calcs
% b4dot: Transitionary term for partial calcs
% b5dot: Transitionary term for partial calcs
% b6dot: Transitionary term for partial calcs
% dxda: Partial derivative of x wrt a
% dxdlam: Partial derivative of x wrt lambda
% dxdh: Partial derivative of x wrt h
% dxdk: Partial derivative of x wrt k
% dxdp: Partial derivative of x wrt p
% dxdq: Partial derivative of x wrt q
% dxdotda: Partial derivative of xdot wrt a
% dxdotdlam: Partial derivative of xdot wrt lambda
% dxdotdh: Partial derivative of xdot wrt h
% dxdotdk: Partial derivative of xdot wrt k
% dxdotdp: Partial derivative of xdot wrt p
% dxdotdq: Partial derivative of xdot wrt q
% 
% OUTPUTS = 
% J: Jacobian between the Equinoctial element set and the PQW frame, depeneding on the direction of choice
% =======================================================================

a = COE(1); e = COE(2); i = COE(3); O = COE(4); w = COE(5); M_t = COE(6);
P = [1; 0; 0]; 
Q = [0; 1; 0]; 
R = [0; 0; 1];

[E,~] = MAtoEccTA(M_t,e);
r = a * (1 - e * cos(E));
n = sqrt(mu/a^3);

X = a * (cos(E) - e); 
Y = a * sqrt(1-e^2) * sin(E);
Xdot = -n * a^2 * sin(E) / r; 
Ydot = n * a^2 * sqrt(1 - e^2) * cos(E) / r;
x = P * X + Q * Y; 
xdot = P * Xdot + Q * Ydot;

L = a^2  * (e * cos(E) - 1 - sin(E)^2) / r; 
M = a^2 * sin(E) * (cos(E) - e) / r / sqrt(1-e^2);
Ldot = n * a^4 * (e - 2 * cos(E) + e * cos(E)^2) * sin(E) / r^3;
Mdot = n * a^4 * (e^2 - 1 - e * cos(E) + 2 * cos(E)^2 - e * cos(E)^3) / r^3 / sqrt(1 - e^2);

b1 = L * P + M * Q; 
b2 = 1/e * (Q * X - P * Y - xdot/n); 
b3 = Q * X - P * Y;
b4 = (X * sin(w + O) + Y * cos(w + O)) * R; 
b5 = (X * cos(w + O) - Y * sin(w + O)) * R;
b6 = (X * sin(w) + Y * cos(w)) * R; 
b1dot = Ldot * P + Mdot * Q;
b2dot = 1/e * (Q * Xdot - P * Ydot + n * (a/r)^3 * x); 
b3dot = Q * Xdot - P * Ydot;
b4dot = (Xdot * sin(w + O) + Ydot * cos(w + O)) * R;
b5dot = (Xdot * cos(w + O) - Ydot * sin(w + O)) * R;
b6dot = (Xdot * sin(w) + Ydot * cos(w)) * R;

dxda  = 1/a * (x - 3/2 * xdot * t); 
dxdlam = xdot / n;
dxdh = b1 * sin(w + O) + b2 * cos(w + O);
dxdk = b1 * cos(w + O) - b2 * sin(w + O);
dxdp = -cos(O) * sin(i) * b3 - (1 + cos(i)) * b5;
dxdq = sin(O) * sin(i) * b3 + (1 + cos(i)) * b4;
dxdotda = -(1/(2*a)) * (xdot - 3 * mu * (x * t / r^3));
dxdotdlam = -n * (a/r)^3 * x;
dxdotdh = b1dot * sin(w + O) + b2dot * cos(w + O);
dxdotdk = b1dot * cos(w + O) - b2dot * sin(w + O);
dxdotdp = -cos(O) * sin(i) * b3dot - (1 + cos(i)) * b5dot;
dxdotdq = sin(O) * sin(i) * b3dot + (1 + cos(i)) * b4dot;


if dir == 1
    J = [dxda dxdh dxdk dxdp dxdq dxdlam;...
        dxdotda dxdotdh dxdotdk dxdotdp dxdotdq dxdotdlam];
else
    J = [dxda dxdh dxdk dxdp dxdq dxdlam;...
        dxdotda dxdotdh dxdotdk dxdotdp dxdotdq dxdotdlam] \ eye(6);
end

end