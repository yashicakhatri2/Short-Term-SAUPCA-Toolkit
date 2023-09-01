function [constants, ICState] = constantsAndInitialState()
% Define constants and initial states of the two objects to set up the
% conjunction.
% 
% IMPORTANT NOTE: Variables marked with * prefix in introductory comments 
% can be modified % by usr without a need for any modification in the code.
% ========================================================================
% VARIABLES = 
% *r1: [km] Hard body radius of object 1
% *r2: [km] Hard body radius of object 2 
% *a: Semi-major axis  
% *e: Eccentricity
% *i: Inclination
% *O: Right Ascension of the Ascending Node
% *w: Arguemnt of periapsis
% *M: Mean Anomaly
% *P: Cartesian covariance for object 1
% *P2: Cartesian covariance for object 2
% Eq1: Equinoctial element set for object 1
% State1: Cartesian state for object 1
% Del1: Delaunay set for object 1
% State2: Cartesian state for object 2
% Del2: Delaunay set for object 2
% ICState: Initial states of the two objects
%   obj1: Parameters related to object 1
%       COE: Classical orbital element set
%       Eq: Equinoctial element set
%       Del: Delaunay set
%       Cart: Cartesian state
%       PoCart: Cartesian covariance
%       PoPoincare: Poincare covariance
%   obj2: Parameters related to object 2 (same sub-parameters as obj1)
% constants:
%   j2: Earth J2 
%   mu: [km^3/s^2] Earth gravitational parameter assuming negligible satellite mass
%   rE: [km] Radius of the Earth
%   MU: [km^3/s^2] Sun gravitational parameter
%   *rho: Reflectivity of the satellite
%   *Aom: [km^2/kg] Area over mass ratio of the satellite
%   Psp: [kg/s^2/km] Solar constant
%   SunT: Period of 1 Earth revolution around the Sun
%   R_star: Combined hard body radii of objects 1 and 2 (r1 + r2)
%   mu_units: Dimensionless Earth gravitational parameter
%   MU_units: Dimensionless Sun gravitational parameter
%   *usePredefinedCases: 1 - Flag to use predefined cases, 0 - defined your own test case
%   *case_flag: Test case number from define_cases function 
%   *points: Number of Monte Carlo points currently being evaluated for each object
%   *JMAX: DOES NOT WORK WITH LESS THAN ~15. List of number of GMM components to split each object into
%   *par: Flag for parallel runs for high # points, 1 - allow parallel runs
%   *nodes: Choose number of parallel nodes
%   Both_GMM: Allow both objects to have a distribution that gets split into GMMs
%   STTOrder: Order of the STT used in the semi-analytical propagation
%   StateDim: Dimension of the state, keep 8 because additional state elements get added for the augmented Hamiltonian
%   options: Option to define the ode propagation tolerances
%   *plot_GMM_ell: 1 - plot GMM ellipses (usually used when testFlagFig2 = 1 and testing at a set time t)
%   *testFig2: Plot - Testing the GMM spread at a set time t = P_prop from function define cases. JMAX FIXED AT 15 COMPONENTS. DOES NOT PROVIDE ACCURATE PROBABILITY OF COLLISION RESULTS.
%   *testFig4: Plot - Check the weight + mean distribution of the GMMs. 'runGMMSTTMethod' = 1 required
%   STTDim: One dimension of the n-dimensional State Transition Tensor
%   initiateSTT: Initiate the STTs with an identity matrix in the first order, and 0s elsewhere
%   *P_prop: Nominal time of closest approach (TCA)
%   JMAX: Number of GMM components chosen for each object in conjunction
% flagNewSymbolicCalcs: **** KEEP AT 0 ***** Flag to recompute STT Equations of Motion * ONLY TOGGLE to 1 IF NEW DYNAMICS ARE INTRODUCED.
% newSTTEOMs: Recomputed Equations of Motion for STTs
% PQW_to_IJK: 3x3 Rotation matrix from IJK to PQW frame for object 1
% PQW_to_IJK2: 3x3 Rotation matrix from IJK to PQW frame for object 1
% J1: 6x6 Rotation matrix from IJK to PQW frame for object 1
% J1_2: 6x6 Rotation matrix from IJK to PQW frame for object 2
% J2: 6x6 Jacobian to transform a matrix from PQW to poincare frame
% J2_2:  6x6 Jacobian to transform a matrix from PQW to poincare frame
% STTo: Initial State Transition Matrix (STM) 
%
% OUTPUTS = 
% constants: Defined above
% ICState: Defined above
% ========================================================================

constants.j2 = 1.08262617385222E-3;
constants.mu = 398600.4415; 
constants.rE = 6378.137;
constants.MU = 1.32712428E11;
constants.rho = 0.2; 
constants.Aom = 2E-6; 
constants.Psrp = 4.57E-3; 
constants.SunT = 365.25 * 24 * 3600; 
r1 = 20 / 1000;
r2 = 20 / 1000;
constants.R_star = r1 + r2;
constants.mu_units = constants.mu / constants.rE^3 * 3600^2;
constants.MU_units = constants.MU / constants.rE^3 * 3600^2;
constants.usePredefinedCases = 1;

flagNewSymbolicCalcs = 0; % ONLY UPDATE IF THERE IS AN UPDATE IN CHOSEN DYNAMICS

% Compute New Jacobian for STT EOM calcs
% ONLY UPDATE IF THERE IS AN UPDATE IN CHOSEN DYNAMICS
if flagNewSymbolicCalcs == 1
    newSTTEOMs = SymbolicSDSJ2Computer(constants.mu, constants.rE, constants.j2, 2, "SP", 0);
end

% Object 1 Specifications
a = 8000; e = 0.15; i = 60 * pi / 180; O = 0; w = 0; M = 0; ICState.obj1.COE = [a;e;i;O;w;M];
Eq1 = COE_to_Equinoctial(ICState.obj1.COE, constants.mu, 1);
State1 = COE_to_Cartesian(ICState.obj1.COE, constants.mu, 1); % km and rad input
Del1 = COE_to_Delaunay(ICState.obj1.COE,constants.mu,1);
ICState.obj1.Eq = Eq1; ICState.obj1.Cart = State1; ICState.obj1.Del = Del1;

% Chose Test Types + Propagation Time
constants.case_flag = 1;
constants.points = 1E2;
constants.JMAX = [15];
constants.par = 0; constants.nodes = 4; 
constants.Both_GMM = 1;
constants.STTOrder = 2;
constants.StateDim = 6;
constants.STTDim = 8;
constants.options = odeset('RelTol', 1E-13,'AbsTol',1E-13);

constants.testFig2 = 0; constants.plot_GMM_ell = 0;
constants.testFig4 = 0; 

% Object 2 orbit specification
if constants.usePredefinedCases == 1
    [ICState.obj2.COE, constants.P_prop] = define_cases(constants.case_flag); % Choses the test case
else
    ICState.obj2.COE = [9843.33683736869         0.241038139941096          1.07980901556839        -0.191768531625329         0.970045507875877         -1.52824344879098]';
    constants.P_prop = 1.5 * 24 * 3600;
end
ICState.obj2.Eq = COE_to_Equinoctial(ICState.obj2.COE, constants.mu, 1);
Del2 = COE_to_Delaunay(ICState.obj2.COE,constants.mu,1);
State2 = Equinoctial_to_Cartesian(ICState.obj2.Eq, constants.mu);
ICState.obj2.Cart = State2; ICState.obj2.Del = Del2;

% Other Constants definition
if constants.testFig2 == 1; constants.JMAX = [15]; end % Replaces given JMAX (#GMM comps) with a small JMAX when testing a small # GMMs at set time t = P_prop

% Initialise the State Transition Tensor
STTo{1} = eye(constants.STTDim);
STTo{2} = zeros(constants.STTDim,constants.STTDim,constants.STTDim);
constants.initiateSTT = STTo;

% New State + Covariance defintions
P = [0.0067664 -0.0029183 0.0027112 -9.9816E-7 -1.7636E-7 2.1797E-6;...
    -0.0029183 0.005348 -0.0011671 -1.5861E-6 -3.5203E-7 3.3414E-6;...
    0.0027112 -0.0011671 0.001087 -3.9883E-7 -7.5945E-8 8.6148E-7;...
    -9.9816E-7 -1.5861E-6 -3.9883E-7 9.4587E-9 -1.1375E-10 1.1511E-9;...
    -1.7636E-7 -3.5203E-7 -7.5945E-8 -1.1375E-10 9.8844E-9 8.5671E-11;...
    2.1797E-6 3.3414E-6 8.6148E-7 1.1511E-9 8.5671E-11 7.2859E-9];
P2 = P; 
ICState.obj1.PoCart = P;
ICState.obj2.PoCart = P2;

% Poincare Uncertainty + Covariance Definition
PQW_to_IJK = PQW_to_IJK_transform(ICState.obj1.COE,0);
PQW_to_IJK2 = PQW_to_IJK_transform(ICState.obj2.COE,0);
J1 = [PQW_to_IJK zeros(3,3); zeros(3,3) PQW_to_IJK]; J1_2 = [PQW_to_IJK2 zeros(3,3); zeros(3,3) PQW_to_IJK2];
J2 = JacobianCalc(ICState.obj1.COE, constants.mu, 0, 0); J2_2 = JacobianCalc(ICState.obj2.COE, constants.mu, 0, 0);
ICState.obj1.PoPoincare = J2 * J1 * P * J1' * J2';
ICState.obj2.PoPoincare = J2_2 * J1_2 * P2 * J1_2' * J2_2';

end