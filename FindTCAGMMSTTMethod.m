function [so, sf, COE1out, COE2out, GMM_mean1i, numSTT1, numSTT2] = FindTCAGMMSTTMethod(mu_old1, mu_old2, ICState, constants)
% This function finds the time of closest approach between the two objects
% in conjunction using a semi-analytical method of propagation.
% 
% =======================================================================
% INPUTS = 
% mu_old1: GMM component mean for object 1 
% mu_old2: GMM component mean for object 2
% ICState: Initial state of objects 
%   obj1: 
%       COE: Initial state of object 1 in COE
%   obj2: 
%       COE: Initial state of object 2 in COE
% constants:
%   P_prop: Nominal Time of Closest Approach
%   mu: Earth gravitational parameter
%   rE: Radius of the Earth
%   STTOrder: Order of the STT used in the semi-analytical propagation
%   testFig2: Plot - Testing the GMM spread at a set time t = P_prop from function define cases. JMAX FIXED AT 15 COMPONENTS.
%   initiateSTT: Initiate the STTs with an identity matrix in the first order, and 0s elsewhere
%
% VARIABLES = 
% P_nom: Nominal Time of Closest Approach
% propTime1: Propagation start time
% propTime2: Propagation end time
% mu: Earth gravitational parameter
% rE: Radius of the Earth
% fig2: Plot - Testing the GMM spread at a set time t = P_prop from function define cases. JMAX FIXED AT 15 COMPONENTS.
% COE1: Initial COE object 1 at the start of ode propagation
% COE2: Initial COE object 2 at the start of ode propagation
% numericSTT1: Initial STT of object 1 at the ode propagation start time
% numericSTT2: Initial STT of object 2 at the ode propagation start time
% iter: Iteration counter
% t: Time update to the actual TCA
% Out1: Propagation output from propagateWithDynamics for object 1
% Out2: Propagation output from propagateWithDynamics for object 2
% Off1: Final de-normalised offset for object 1
% Off2: Final de-normalised offset for object 2
% Del_updatedi1: Object 1 mean state at the final time 
% Del_updatedi2: Object 2 mean state at the final time 
% mu_updated1: Mapped GMM component delta mean for object 1
% mu_updated2: Mapped GMM component delta mean for object 2
% mu_new1: De-normalised mu_updated1
% mu_new2: De-normalised mu_updated2
% New_plot1: Delta mean (mu_updated1) added to Object 1 mean state at the end, converted to osculating state 
% New_plot2: Delta mean (mu_updated2) added to Object 2 mean state at the end, converted to osculating state 
% NewOff1: Final offset 1 to convert mean state of the deviated mean to osculating at the final time
% NewOff2: Final offset 2 to convert mean state of the deviated mean to osculating at the final time
% GMM_mean1i: Final GMM component mean for object 1 in Cartesian state
% GMM_mean2i: Final GMM component mean for object 2 in Cartesian state
% pR: Relative position between the two objects at the final time
% vR: Relative velocity between the two objects at the final time
% 
% OUTPUTS = 
% so: Relative position between the two object at the nominal TCA
% sf: Relative position between the two objects at the actual TCA
% COE1out: Final object 1 state at actual TCA in COE
% COE2out: Final object 2 state at actual TCA in COE
% GMM_mean1i: GMM component mean for object 1 in Cartesian state at final time, used for plotting
% numSTT1: Final STTs for object 1
% numSTT2: Final STTs for object 2
% =======================================================================

% Constants
P_nom = constants.P_prop;
propTime1 = 0;
propTime2 = P_nom;
mu = constants.mu;
rE = constants.rE;
fig2 = constants.testFig2;

% Initial Conditions for TCA iterator
COE1 = ICState.obj1.COE;
COE2 = ICState.obj2.COE;
numericSTT1 = constants.initiateSTT;
numericSTT2 = constants.initiateSTT;
iter = 0;

t = 10;
% TCA iterator
while abs(t) > 1E-10 && iter < 50
    iter = iter + 1;

    % STT Calculations at current TCA
    Out1 = propagateWithDynamics(COE1, numericSTT1, constants, [propTime1/3600, propTime2/3600], 1);
    Out2 = propagateWithDynamics(COE2, numericSTT2, constants, [propTime1/3600, propTime2/3600], 1);
    numericSTT1 = Out1{2};
    numericSTT2 = Out2{2};
    numericSTT1{3} = Out1{3};
    numericSTT2{3} = Out2{3};
    Off1 = normalize(Out1{3},constants.rE,'vec','Del',0)';
    Off2 = normalize(Out2{3},constants.rE,'vec','Del',0)';
    Del_updatedi1 =  Out1{1}(1:6)' - Off1(1:6);
    Del_updatedi2 =  Out2{1}(1:6)' - Off2(1:6);

    % Propagating object 1
    mu_updated1 = STTcalcs2BP(mu_old1, constants.STTOrder, constants, numericSTT1);
    mu_new1(:,1) = normalize(mu_updated1,rE,'vec','Del',0);
    COE1 = COE_to_Delaunay(Out1{1}(1:6)',mu,0);
    New_plot1 = Del_updatedi1 + mu_new1;
    NewOff1 = getOffset([New_plot1;0;0],propTime2,constants,1);
    New_plot1 = New_plot1 + NewOff1(1:6);
    GMM_mean1i = COE_to_Cartesian(COE_to_Delaunay(New_plot1,mu,0),mu,1);
    
    % Propagating object 2 
    mu_updated2 = STTcalcs2BP(mu_old2, constants.STTOrder, constants, numericSTT2);
    mu_new2(:,1) = normalize(mu_updated2,rE,'vec','Del',0);
    COE2 = COE_to_Delaunay(Out2{1}(1:6)',mu,0);
    New_plot2 = Del_updatedi2 + mu_new2;
    NewOff2 = getOffset([New_plot2;0;0],propTime2,constants,1);
    New_plot2 = New_plot2 + NewOff2(1:6);
    GMM_mean2i = COE_to_Cartesian(COE_to_Delaunay(New_plot2,mu,0),mu,1);
    
    if iter == 1
        so = GMM_mean2i - GMM_mean1i;
    end
    if fig2 == 1
        break
    end

    pR = GMM_mean1i(1:3) - GMM_mean2i(1:3);
    vR = GMM_mean1i(4:6) - GMM_mean2i(4:6);

    t = - dot(pR, vR) / norm(vR)^2;
    propTime1 = propTime2;
    propTime2 = propTime2 + t;
    
end
sf = GMM_mean2i - GMM_mean1i;
COE1out = COE_to_Delaunay(New_plot1, mu, 0);
COE2out = COE_to_Delaunay(New_plot2, mu, 0);
numSTT1 = numericSTT1(1:2);
numSTT2 =  numericSTT2(1:2);
if iter > 45
    if abs(t) > 1E-1
        return
    end
end
end
