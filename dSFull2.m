function dSFull2 = EOM(t, S, constants)
% Full dynamics propagation equations
% =======================================================================
% INPUTS = 
% t: Current time
% S: Current state [position;velocity]
% constants:
%   MU: Gravitational parameter of the Sun
%   j2: Gravitational potential constant J2
%   mu: Gravitational parameter of the Earth
%   Aom: Area over mass ratio of sattelite
%   Psrp: Solar constant
%   rE: Radius of the Earth
%   rho: Reflectivity of the satellite
%
% VARIABLES:
% rE: Radius of the Earth
% J2: Gravitational potential constant J2  
% mu: Gravitational parameter of the Earth
% r: Norm of position vector
% rvec: Position vector
% rI: Position vector component 1
% rJ: Position vector component 2
% rK: Position vector component 3
% a2BP: Acceleration due to 2 body problem
% SSEvec: Sun state wrt Earth
% rrelvec: Relative state between the Sun and the satellite
% rrel: Norm of rrelvec
% b: Solar perturbation strength
% aSRP: Acceleration due to Solar Radiation Pressure
% aI: Acceleration due to J2 gravitational harmonics component 1
% aJ: Acceleration due to J2 gravitational harmonics component 2
% aK: Acceleration due to J2 gravitational harmonics component 3
% aJ2: Acceleration due to J2 gravitational harmonics
% atot: Acceleration due to all perturbative forces combined
% 
% OUTPUT =
% dSFull2: Final velocity + acceleration = ode propagation vector
% =======================================================================

rvec = S(1:3); r = norm(rvec);
rI = rvec(1); rJ = rvec(2); rK = rvec(3);
mu = constants.mu;
rE = constants.rE;
J2 = constants.j2;

% 2BP
a2BP = - mu / r^3 * rvec;

% SRP
SSEvec = COE_to_Cartesian_2(getSun(constants,t),constants.MU,1);
rrelvec = SSEvec(1:3) - rvec; rrel = norm(rrelvec);
b = (1+constants.rho) * constants.Psrp * constants.Aom;
aSRP = -b * rrelvec / rrel;

% J2
aI = -3 * J2 * mu * rE^2 * rI / 2 / r^5 * (1 - 5 * rK^2 / r^2);
aJ = -3 * J2 * mu * rE^2 * rJ / 2 / r^5 * (1 - 5 * rK^2 / r^2);
aK = -3 * J2 * mu * rE^2 * rK / 2 / r^5 * (3 - 5 * rK^2 / r^2);
aJ2 = [aI;aJ;aK];

% Final acceleration
atot = a2BP + aSRP + aJ2;

dSFull2(1:3,1) = S(4:6);
dSFull2(4:6,1) = atot;

end