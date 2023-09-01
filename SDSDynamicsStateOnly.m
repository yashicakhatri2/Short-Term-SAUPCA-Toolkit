function dS = SDSDynamicsStateOnly(t,S,constants)
% SDS Propagation Method for state only (not STTs).
%
% ========================================================================
% SUBFUNCTIONS GENERATED USING MATLAB SYMBOLIC TOOLBOX:
% MeanDynamicsFunction: Coded equations for mean dynamics generated using a Simplified Dynamics System (SDS), coded using MATLAB symbolic toolbox using script: SDS Generator
% 
% =======================================================================
% INPUTS = 
% t: Time
% S: State only
% constants:
%   MU: Gravitational parameter of the Sun
%   MU_units: Normalised gravitational parameter of the Sun
%   j2: Gravitational potential constant J2
%   mu_units: Normalised gravitational parameter of the Earth
%   Aom: Area over mass ratio of sattelite
%   Psrp: Solar constant
%   rE: Radius of the Earth
%   rho: Reflectivity of the satellite
% 
% VARIABLES = 
% l: Delaunay set element 1
% g: Delaunay set element 2
% h: Delaunay set element 3
% L: Delaunay set element 4
% G: Delaunay set element 5
% H: Delaunay set element 6
% ks: New Delaunay set element 7 generalised coordinate for augmented Hamiltonian
% Ks: New Delaunay set element 8 conjugtate momentum for augmented Hamiltonian
% SunCOE: Classical orbital elements of the Sun at this time
% SunDel: Delaunay set of the Sun at this time
% hs: Sun Delaunay set element 3
% cs: cos(Inclination of the Sun)
% sSun: sin(Inclination of the Sun)
% nu: Mean motion of the sun
% RE: Normalised radius of the Earth = 1
% J2: Gravitational potential constant J2  
% mu: Normalised gravitational parameter of the Earth
% Aom: Normalised Area over mass ratio of sattelite 
% Psrp: Normalised Solar constant
% rho: Reflectivity of the satellite
% b: Solar perturbation strength
% 
% OUTPUT =
% dSSTTpropagation: Equations of motion vector for the state ode propagation
% =======================================================================

l = wrapToPi(S(1)); g = wrapToPi(S(2)); h = wrapToPi(S(3)); L = S(4); G = S(5); H = S(6); ks = S(7); Ks = S(8);

SunCOE = getSun(constants,t*3600);
SunDel = normalize(COE_to_Delaunay(SunCOE,constants.MU,1),constants.rE,'vec','Del',1); % Using normalized mu for sun propagation, we only care about h from the Sun coordinates.
hs = SunDel(3); % rad
cs = SunDel(6) / SunDel(5);
sSun = sqrt(1 - cs^2);
nu = constants.MU_units^2 / SunDel(4)^3;

% Using normalized constants for propagation
RE = 1;
J2 = constants.j2;
mu = constants.mu_units;
Aom = constants.Aom / constants.rE^2; % km^2/kg 
Psrp = constants.Psrp * constants.rE * 3600^2; % kg/s^2/km
rho = constants.rho;
b = (1+rho) * Aom * Psrp;

dS = MeanDynamicsFunction(G,H,J2,L,RE,b,cs,g,h,hs,ks,mu,nu,sSun);
end
