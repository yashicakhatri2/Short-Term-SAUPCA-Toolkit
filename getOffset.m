function OUT = getOffset(S, t, constants, flagIF)
% This function computes the offset that allows conversion of a state from
% mean to osculating and vice-versa.
%
% ========================================================================
% SUBFUNCTIONS GENERATED USING MATLAB SYMBOLIC TOOLBOX:
% getInitialDelOffset: Initial offset to convert from osculating state to mean state
% getFinalDelOffset: Final offset to convert from mean state to osculating state
%
% ========================================================================
% INPUTS = 
% t: Current time at which offset needs to be computed
% S: Delaunay state at time t
% constants: 
%   MU: Gravitational parameter of the Sun
%   j2: Gravitational potential constant J2
%   mu_units: Normalised gravitational parameter of the Earth
%   Aom: Area over mass ratio of sattelite
%   Psrp: Solar constant
%   rE: Radius of the Earth
%   rho: Reflectivity of the satellite
%   STTDim: One dimension of the n-dimensional State Transition Tensor
%   options: Option to define the ode propagation tolerances
%   mu: Earth gravitational parameter
%   rE: Radius of the Earth
%   STTOrder: Order of the STT
% flagIF: 0 - initial offset from osculating to mean, 1 - final offset from osculating to mean
%  
% VARIABLES = 
% E: Arbitrary small value
% S: Normalised Delaunay input state
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
% a: Semi-major axis
% e: Eccentricity
% RE: Normalised radius of the Earth = 1
% J2: Gravitational potential constant J2  
% mu: Normalised gravitational parameter of the Earth
% Aom: Normalised Area over mass ratio of sattelite 
% Psrp: Normalised Solar constant
% rho: Reflectivity of the satellite
% b: Solar perturbation strength
% Ecc: Eccentric anomaly
% TA: True anomaly
% r: Distance from Earth
% SOUT: Normalised offset at the current time in the chosen direction
% 
% OUTPUTS = 
% OUT: De-normalised offset Delaunay elements
% ========================================================================

E = 1;
S = normalize(S,constants.rE,'vec','Del',1);
l = wrapToPi(S(1)); g = wrapToPi(S(2)); h = wrapToPi(S(3)); L = S(4); G = S(5); H = S(6); ks = S(7); Ks = S(8);

SunCOE = getSun(constants,t*3600);
SunDel = normalize(COE_to_Delaunay(SunCOE,constants.MU,1),constants.rE,'vec','Del',1); % Only h needed from the Sun coordinates.
hs = SunDel(3); % radians
cs = SunDel(6) / SunDel(5);
sSun = sqrt(1 - cs^2);

% Using normalized constants for propagation
RE = 1;
J2 = constants.j2;
mu = constants.mu_units;
Aom = constants.Aom / constants.rE^2;
Psrp = constants.Psrp * constants.rE * 3600^2; 
rho = constants.rho;

b = (1+rho) * Aom * Psrp;

a = L^2/mu;
e = sqrt(1-(G/L)^2);
[Ecc, TA] = MAtoEccTA(l,e);
r = a * (1 - e * cos(Ecc));
if flagIF == 0
    SOUT = -getInitialDelOffset(E,Ecc,G,H,J2,L,RE,TA,b,cs,g,h,hs,ks,l,mu,r,sSun);
elseif flagIF == 1
    SOUT = -getFinalDelOffset(E,Ecc,G,H,J2,L,RE,TA,b,cs,g,h,hs,ks,l,mu,r,sSun);
end

OUT = normalize(SOUT,constants.rE,'vec','Del',0);
end