function [OUT1, OUT2] = FindTCAMCMethod(COE1, COE2, P_nom, constants)
% This function computes the actual time of closest (TCA) approach using a while
% loop and iterations that use nominal TCA as the initial guess and update
% the TCA with each iteration depending on a norm function that uses the
% relative state between the two objects. The states at the final TCA are
% the outputs of this function.
%
% ========================================================================
% INPUTS = 
% COE1: Classical orbital element set for object 1
% COE2: Classical orbital element set for object 2 
% P_nom: Nominal time of closest approach (TCA)
% constants: Variable with initially defined constants passed through
%   mu: Earth gravitational constant
%   fig2: Flag to plot propagated trajectories at the nominal TCA (skips the TCA finder and outputs values at nominal TCA)
%   
% VARIABLES = 
% mu: Earth gravitational constant
% t: Actual TCA finder variable gets updated with each iteration of the while loop
% iter: Iteration counter for the while loop
% propTime1: Propagation start time
% propTime2: Propagation end time
% Offset1: Final offset for object 1 at the current TCA in the while loop
% Offset2: Final offset for object 2 at the current TCA in the while loop
% COE1: Classical orbital element set for object 1 at the final propagation time at each while loop iteraiton
% COE2: Classical orbital element set for object 2 at the final propagation time at each while loop iteraiton
% S1: Cartesian state for object 1 at the final propagation time at each while loop iteraiton
% S2: Cartesian state for object 1 at the final propagation time at each while loop iteraiton
% sR: Relative cartesian state between objects 1 and 2 (S1 - S2)
%
% OUTPUT = 
% OUT1: Cartesian state for object 1 at the final calculated TCA
% OUT2: Cartesian state for object 2 at the final calculated TCA
% ========================================================================

mu = constants.mu; 
t = 10;
iter = 0;
propTime1 = 0;
propTime2 = P_nom;
Offset1 = 0;
Offset2 = 0;

while abs(t) > 1E-10 && iter < 20
    clearvars S1 S2
    iter = iter + 1;

    Del1 = propagateWithDynamics(COE1, Offset1, constants, [propTime1/3600, propTime2/3600], 0);
    Del2 = propagateWithDynamics(COE2, Offset2, constants, [propTime1/3600, propTime2/3600],  0);
    Offset1 = Del1{2};
    Offset2 = Del2{2}; 
    COE1 = COE_to_Delaunay(Del1{1},mu,0);
    COE2 = COE_to_Delaunay(Del2{1},mu,0);
    S1 = COE_to_Cartesian(COE1,mu,1);
    S2 = COE_to_Cartesian(COE2,mu,1);
    
    sR = S1 - S2;
    t = - dot(sR(1:3), sR(4:6)) / norm(sR(4:6))^2;
    propTime1 = propTime2;
    propTime2 = propTime2 + t;

    if constants.testFig2 == 1
        break
    end
end

OUT1 = S1;
OUT2 = S2;

if iter > 17
    if abs(t) > 1
        return
    end
end

end