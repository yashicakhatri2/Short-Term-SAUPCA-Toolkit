function [PC_GMMi, GMM_MEAN, Q_new] = ComputeGMMSTTMethodPc(ICState,constants)
% Monte Carlo Collision of Probability of Collision calculations.
% =======================================================================
% INPUTS = 
% ICState: Initial states of the two objects
%   obj1: Parameters related to object 1
%       Cart: Cartesian state
%       PoCart: Cartesian covariance
%   obj2: Parameters related to object 2 (same sub-parameters as obj1)
% constants:
%   mu: [km^3/s^2] Earth gravitational parameter assuming negligible satellite mass
%   R_star: Combined hard body radii of objects 1 and 2 (r1 + r2)
%   points: Number of Monte Carlo points currently being evaluated for each object
%   par: Flag for parallel runs for high # points, 1 - allow parallel runs
%   nodes: Choose number of parallel nodes
%   Both_GMM: Allow both objects to have a distribution that gets split into GMMs
%   testFig2: Plot - Testing the GMM spread at a set time t = P_prop from function define cases. JMAX FIXED AT 15 COMPONENTS.
%   P_prop: Nominal time of closest approach (TCA)
%
% VARIABLES = 
% DelNew1: Delaunay state of  object 1 extended to include 2 additional elements for the Augmentation of the Hamiltonian, with offset added for mean state
% DelNew2: Delaunay state of  object 2 extended to include 2 additional elements for the Augmentation of the Hamiltonian, with offset added for mean state
% corrState: The corrected (mean) states defined for both objects for use in the semi-analytical propagation
%   obj1: Parameters related to object 1
%       Del: Mean Delaunay state     
%       COE: Mean Classical Orbital element set
%       Eq: Mean Equinoctial element set
%       PoCart: Cartesian covariance with corrected (mean) state used for transformation
%   obj2: Parameters related to object 2 (same sub-parameters as obj1)
% PQW_to_IJK: 3x3 Rotation matrix from IJK to PQW frame for object 1
% PQW_to_IJK2: 3x3 Rotation matrix from IJK to PQW frame for object 1
% J1: 6x6 Rotation matrix from IJK to PQW frame for object 1
% J1_2: 6x6 Rotation matrix from IJK to PQW frame for object 2
% J2: 6x6 Jacobian to transform a matrix from PQW to poincare frame
% J2_2:  6x6 Jacobian to transform a matrix from PQW to poincare frame
% jmax1: Number of GMM components for object 1
% jmax2: Number of GMM components for object 2
% mean_new1: GMM component means for object 1
% Q_bar1: GMM component covariance for object 1
% w_bar1: GMM component weights for object 1
% mean_new2: GMM component means for object 2
% Q_bar2: GMM component covariance for object 2
% w_bar2: GMM component weights for object 2
% COE_mean_1: mean_new1 converted to Classical Orbital Elements
% COE_mean_2: mean_new2 converted to Classical Orbital Elements
% mean1_del: COE_mean_1 coverted to Delaunay
% mean2_del: COE_mean_2 coverted to Delaunay
% mu_old1i_inter: Delta Mean object 1 = Individual Component Mean - Nominal Trajectory
% mu_old2i_inter: Delta Mean object 1 = Individual Component Mean - Nominal Trajectory
% mu_old1i: Normalised mu_old1i_inter = Main component mean for object 1 to be propagated
% mu_old2i: Normalised mu_old2i_inter = Main component mean for object 2 to be propagated
% Eq_to_Del1: Rotation Jacobian to convert object 1 component covariance to Dealunay Frame
% Eq_to_Del2: Rotation Jacobian to convert object 2 component covariance to Dealunay Frame
% Q_bar1_inter: Delaunay object 1 component covariance
% Q_bar2_inter: Delaunay object 2 component covariance
% Q_bar1_norm: Normalised Delaunay covariance for object 1 GMM components - to be propagated
% Q_bar2_norm: Normalised Delaunay covariance for object 2 GMM components - to be propagated
% PC_GMM: Total semi-analyticla probability of collision
% mu_old1: GMM component mean for current object 1 in the nested for loops
% mu_old2: GMM component mean for current object 2 in the nested for loops
% so: Difference in object states at nominal Time of Closest Approach in encounter plane
% sf: Difference in object states at actual Time of Closest Approach in encounter plane
% New_plot_COE: Final COE set for object 1 at actual TCA
% New_plot_COE_2: Final COE set for object 2 at actual TCA
% GMM_mean: GMM means propagated to actual TCA for plotting
% numericSTT1: STTs for object 1 at actual TCA used for covariance propagation
% numericSTT2: STTs for object 2 at actual TCA used for covariance propagation
% Q_new_1: Propagated object 1 component covariance
% Q_new_2: Propagated object 2 component covariance
% PC_GMMith: Probability of collision for ith and jth GMM componented of the objects 1 and 2 respectively
% 
% OUTPUT =
% PC_GMMi: Cummulative probability of collision for qth # of GMM components
% GMM_MEAN: GMM means propagated to actual TCA for plotting
% Q_new: Propagated object 1 component covariance for plotting
% =======================================================================

DelNew1 = [ICState.obj1.Del; 0; 0] + getOffset([ICState.obj1.Del;0;0],0,constants,0);
corrState.obj1.Del = DelNew1(1:6);
corrState.obj1.COE = COE_to_Delaunay(corrState.obj1.Del,constants.mu,0);
corrState.obj1.Eq = COE_to_Equinoctial(corrState.obj1.COE,constants.mu,1);

DelNew2 = [ICState.obj2.Del; 0; 0] + getOffset([ICState.obj2.Del;0;0],0,constants,0);
corrState.obj2.Del = DelNew2(1:6);
corrState.obj2.COE = COE_to_Delaunay(corrState.obj2.Del,constants.mu,0);
corrState.obj2.Eq = COE_to_Equinoctial(corrState.obj2.COE,constants.mu,1);

% Poincare Uncertainty + Covariance Definition
PQW_to_IJK = PQW_to_IJK_transform(corrState.obj1.COE,0);
PQW_to_IJK2 = PQW_to_IJK_transform(corrState.obj2.COE,0);
J1 = [PQW_to_IJK zeros(3,3); zeros(3,3) PQW_to_IJK]; J1_2 = [PQW_to_IJK2 zeros(3,3); zeros(3,3) PQW_to_IJK2];
J2 = JacobianCalc(corrState.obj1.COE, constants.mu, 0, 0); J2_2 = JacobianCalc(corrState.obj2.COE, constants.mu, 0, 0);
corrState.obj1.PoPoincare = J2 * J1 * ICState.obj1.PoCart * J1' * J2';
corrState.obj2.PoPoincare = J2_2 * J1_2 * ICState.obj2.PoCart * J1_2' * J2_2';

%% Loop for GMM-STT Method   
for q = 1:length(constants.JMAX)
    % Set current loop to run for current number of GMM components
    jmax1 = constants.JMAX(q);
    [mean_new1, Q_bar1, w_bar1] = GMM_Component_Calcs(jmax1, corrState.obj1.Eq, corrState.obj1.PoPoincare,constants.testFig4); % Gives the GMM distribution means, covs and weights
    jmax2 = jmax1;
    [mean_new2, Q_bar2, w_bar2] = GMM_Component_Calcs(jmax2, corrState.obj2.Eq, corrState.obj2.PoPoincare,constants.testFig4); % Gives the GMM distribution means, covs and weights
    
    %% Convert mean and covariance to Delaunay
    mu_old1i = NaN(6,jmax1);
    mu_old2j = NaN(6,jmax2);
    for i = 1:size(mean_new1,2)
        COE_mean_1 = COE_to_Equinoctial(mean_new1(:,i), constants.mu, 0);
        mean1_del = [COE_to_Delaunay(COE_mean_1,constants.mu,1); 0; 0];
        COE_mean_2 = COE_to_Equinoctial(mean_new2(:,i), constants.mu, 0);
        mean2_del = [COE_to_Delaunay(COE_mean_2,constants.mu,1); 0; 0];
        mu_old1i_inter = mean1_del(1:6) - corrState.obj1.Del;
        mu_old2j_inter = mean2_del(1:6) - corrState.obj2.Del;
        mu_old1i(:,i) = normalize(mu_old1i_inter,constants.rE,'vec','Del',1); 
        mu_old2j(:,i) = normalize(mu_old2j_inter,constants.rE,'vec','Del',1); 
    end
    % Covariance
    Eq_to_Del1 = computeEqDelPartials(constants, corrState.obj1.COE, 1);
    Eq_to_Del2 = computeEqDelPartials(constants, corrState.obj2.COE, 1);
    Q_bar1_inter = Eq_to_Del1 * Q_bar1 * Eq_to_Del1';
    Q_bar2_inter = Eq_to_Del2 * Q_bar2 * Eq_to_Del2';
    Q_bar1_norm = normalize(Q_bar1_inter,constants.rE,'mat','Del',1); 
    Q_bar2_norm = normalize(Q_bar2_inter,constants.rE,'mat','Del',1); 
    
 
    %% GMM Pc Calculations
    PC_GMM = 0;
    GMM_MEAN = zeros(6,jmax1);
    if constants.par == 0 
        for i = 1:jmax1
            for j = 1:jmax2
                clearvars mu_old1 mu_old2 so sf New_plot_COE New_plot_COE_2 P_prop_2 GMM_mean numericSTT1 numericSTT2 Q_new_1 Q_new_2
                mu_old1 = mu_old1i(:,i);            
                mu_old2 = mu_old2j(:,j);
                [so, sf, New_plot_COE, New_plot_COE_2, GMM_mean, numericSTT1, numericSTT2] = FindTCAGMMSTTMethod(mu_old1, mu_old2, ICState, constants);

                [Q_new_1, Q_new_2] = Covariance_Propagation(Q_bar1_norm, Q_bar2_norm, New_plot_COE, New_plot_COE_2, mu_old1, mu_old2, constants, numericSTT1, numericSTT2); % Propagates both covariances

                PC_GMMith = GMM_Pc_Cacls_Fun(so, sf, Q_new_1+Q_new_2, constants.R_star); % GMM Probability of Collision Calculator
                if constants.testFig2 == 1; GMM_MEAN(:,i) = GMM_mean; else GMM_mean = 1; end
                if (constants.testFig2 == 1 && constants.plot_GMM_ell == 1); Q_new{i} = Q_new_1; else Q_new = 1; end
                PC_GMM = PC_GMM + (w_bar1(i) * w_bar2(j)) * PC_GMMith; % Total Pc calc: Sums up the weighed Pc for each individual GMM component
            end
        end
    else
        parpool('local',constants.nodes); 
        parfor i = 1:jmax1
            for j = 1:jmax2
                mu_old1 = mu_old1i(:,i);            
                mu_old2 = mu_old2j(:,j);
                [so, sf, New_plot_COE, New_plot_COE_2, ~, numericSTT1, numericSTT2] = FindTCAGMMSTTMethod(mu_old1, mu_old2, ICState, constants); 
                
                [Q_new_1, Q_new_2] = Covariance_Propagation(Q_bar1_norm, Q_bar2_norm, New_plot_COE, New_plot_COE_2, mu_old1, mu_old2, constants, numericSTT1, numericSTT2); % Propagates both covariances

                PC_GMMith = GMM_Pc_Cacls_Fun(so, sf, Q_new_1+Q_new_2, constants.R_star); % GMM Probability of Collision Calculator
                PC_GMM = PC_GMM + (w_bar1(i) * w_bar2(j)) * PC_GMMith; % Total Pc calc: Sums up the weighed Pc for each individual GMM component
            end
        end
        GMM_mean = 1; Q_new = 1;
        delete(gcp('nocreate')) 
    end
   
    PC_GMMi(q) = PC_GMM; % Places total Pc according to the GMM # comps vector
end

end
