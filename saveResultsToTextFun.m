function saveResultsToTextFun(ICState, results, constants)
% This function saves the function results to a text file.
%
% =======================================================================
% INPUTS = 
% ICState: Initial states of the two objects
%   obj1: Parameters related to object 1
%       COE: Classical orbital element set
%       PoCart: Cartesian covariance
%   obj2: Parameters related to object 2 (same sub-parameters as obj1)
% results:
%   MC_Pc: Monte Carlo probability of collision
%   points: Number of points used to compute MC_Pc
%   MCRunTime: Runtime for MC_Pc calculations, some are saved, use tic-toc to compute and save new ones
%   GMMSTT_Pc: Semi-analytical probability of collision
%   JMAX: List of number of GMMs components for each object used to compute GMMSTT_Pc
%   GMMSTTRunTime: List of times to compute GMMSTT_Pc using # components listed in JMAX
% constants:
%   case_flag: Test number from define_cases function 
% 
% VARIABLES = 
% MCRunTime: Runtime of the Monte Carlo calculations
% GMMSTTRunTime: Runtime of the semi-analytial method calculations
% fileName: Name of the file to store the results in 
% fileID: Open and write to file
% 
% OUTPUTS = 
% Saves results to a text file.
% =======================================================================

MCRunTime = results.MCRunTime;
GMMSTTRunTime = results.GMMSTTRunTime;

% Print results to file
fprintf('Using dynamics: %s\n','SDS+SP')
fprintf('MC Pc = %.10f \n', results.MC_Pc);
for q = 1:length(results.JMAX)
    fprintf('GMM-STT Method Pc for %i GMMs = %.10f \n', results.JMAX(q), results.GMMSTT_Pc(q) );
end

% Creates a new txt file and saves results
fileName = strcat('testResults_Case',num2str(constants.case_flag),'_Using_','SDS+SP','.txt');
fileID = fopen(fileName,'w');
fprintf(fileID,'Nominal TCA = %.2f days \n\n', constants.P_prop / 3600 / 24);
fprintf(fileID,'Object 1 Specs:\nCOE = [%.2f km, %.2f, %.2f rad, %.2f rad, %.2f rad, %.2f rad] \nInitial Cartesian Covariance [km;km/s] = \n', ICState.obj1.COE); % [a;e;i;O;w;M]
fprintf(fileID,'%.6e, %.6e, %.6e, %.6e, %.6e, %.6e\n', ICState.obj1.PoCart); % [a;e;i;O;w;M]
fprintf(fileID,'\n\nObject 2 Specs:\nCOE = [%.2f km, %.2f, %.2f rad, %.2f rad, %.2f rad, %.2f rad] \nInitial Cartesian Covariance [km;km/s] = \n', ICState.obj2.COE); % [a;e;i;O;w;M]
fprintf(fileID,'%.6e, %.6e, %.6e, %.6e, %.6e, %.6e\n', ICState.obj2.PoCart); % [a;e;i;O;w;M]
fprintf(fileID,'\n\nResults from analysis:\nMC Pc = %.10f for %.2e points, , with runtime of %.2f seconds.\n', results.MC_Pc,results.points, MCRunTime);
for iter = 1:length(results.JMAX)
    fprintf(fileID,'GMM-STT Method Pc for %i GMMs = %.10f, with runtime of %.2f seconds.\n', results.JMAX(iter), results.GMMSTT_Pc(iter), GMMSTTRunTime(iter));
end
fclose(fileID);
end