function OUT = getSun(constants,t)
% Sun geometric propagator.
% ========================================================================
% INPUTS = 
% t: Current time to grab Sun state
% constants: 
%   MU: Gravitational parameter of the Sun
%   SunT: Period of 1 Earth revolution around the Sun
%  
% VARIABLES = 
% T: Period of 1 Earth revolution around the Sun
% a: Semi-major axis of a circular Sun orbit around the Earth
% n: Mean motion of a circular Sun orbit around the Earth
% M: Anomaly at t
%
% OUTPUTS = 
% OUT: Sun COE state at t
% ========================================================================

T = constants.SunT;
a = ((T / 2 / pi)^2 * constants.MU)^(1/3);
n = sqrt(constants.MU/a^3);

M = wrapToPi(n * t);

OUT = [a;0;0;0;0;M];
end
