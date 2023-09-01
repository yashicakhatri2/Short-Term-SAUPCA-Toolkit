% Created and uploaded by: Yashica Khatri August 2023
%
% =======================================================================
% Short-Term Conjunction Analysis Toolkit using a combination of Gaussian
% Mixture Models (GMM) and State Transition Tensors (STTs). The dynamics
% used in this work include J2 and SRP dynamics. Hamiltonian averaging is
% used to achieve the mean dynamics using the chosen perturbation, this is
% called the Simplified Dynamical System (SDS), previoulsy defined by Park
% and Scheeres.
%
% IMPORTANT NOTE: Variables marked with * prefix in introductory comments 
% can be modified % by usr without a need for any modification in the code.
% =======================================================================
% MATLAB toolboxes used:
% optimization_toolbox
% statistics_toolbox
% 
% =======================================================================
% INPUTS = 
% *runMonteCarlo: Flag to run the Monte Carlo simulation
% *runGMMSTTMethod: Flag to run the semi-analytical probability of collision method 
% *plotThings: Flag to plot results
% *saveResultsToText: Save results to an output file
%
% VARIABLES = 
% constants: Initial defined constants
%   points: Number of points initially defined to compute MC_Pc  
%   case_flag: Test number from define_cases function, using predefined test cases
%   JMAX: List of predefined number of GMMs components for each object used to compute GMMSTT_Pc
% results:
%   MC_Pc: Monte Carlo probability of collision
%   points: Number of points used to compute MC_Pc
%   MCRunTime: Runtime for MC_Pc calculations, some are saved, use tic-toc to compute and save new ones
%   GMMSTT_Pc: Semi-analytical probability of collision
%   JMAX: List of number of GMMs components for each object used to compute GMMSTT_Pc
%   GMMSTTRunTime: List of times to compute GMMSTT_Pc using # components listed in JMAX
% ICState: Initial conditions for objects
% GMM_MEAN: GMM component propagated means  
% Q_new: GMM component propagated covariance
% =======================================================================

clc;
clear;
close all;
format longg;
rng('shuffle'); % Shuffle range for random point generation 

%% Initial Setup
% Run Flags
runMonteCarlo = 0; % Run New MC Pc calcs? 
runGMMSTTMethod = 0;  % Run New GMMSTT Method Pc calcs?
plotThings = 1; % Plot results using new or presaved results
saveResultsToText = 0; % Save results to a text file

% Set constants and intitial states of the objects in conjunction
[constants, ICState] = constantsAndInitialState();

%% Probability of Collision Calculations
% MC Pc Calculations
if runMonteCarlo == 1
    results.MC_Pc = Compute_MC_Pc_Cart_1on1(ICState, constants);
    results.points = constants.points;
else
    [results.MC_Pc, results.points, results.MCRunTime] = savedResults(constants.case_flag,1);
end

% GMM-STT Method Pc Calculations
if runGMMSTTMethod == 1
    [results.GMMSTT_Pc, GMM_MEAN, Q_new] = ComputeGMMSTTMethodPc(ICState,constants);
    results.JMAX = constants.JMAX;
else
    [results.GMMSTT_Pc, results.JMAX, results.GMMSTTRunTime] = savedResults(constants.case_flag,2);
    GMM_MEAN = 0; % Test case, have to run full method to plot
    Q_new = 0; % Test case, have to run full method to plot
end

%% Post-processing 
% Plot the Pc comparisons and distribution comparison at t = P_prop
if plotThings == 1
    Main_Plotting(results,constants,GMM_MEAN, Q_new); 
end

% Write results to a text file
if saveResultsToText == 1
    saveResultsToTextFun(ICState, results, constants);
end
