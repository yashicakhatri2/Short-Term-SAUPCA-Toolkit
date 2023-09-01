function [Q_new, Q_new_2] = Covariance_Propagation(Q_bar, Q_bar2, New_plot_COE, New_plot_COE_2,mu_old1,mu_old2, constants, numericSTT1, numericSTT2)
% This function uses STTs to propagate the covariance of both objects to 
% the desired time.
% 
% =======================================================================
% INPUTS = 
% Q_bar: Epoch GMM component covariance for object 1
% Q_bar2: Epoch GMM component covariance covariance for object 2 
% New_plot_COE: Object 1 state in COE at the final propagation time
% New_plot_COE_2: Object 2 state in COE at the final propagation time
% mu_old1: Object 1 GMM component delta mean at epoch 
% mu_old2: Object 2 GMM component delta mean at epoch 
% constants:
%   mu: Earth gravitational parameter
%   rE: Radius of the Earth
%   STTOrder: Order of the STT used in the semi-analytical propagation
% numericSTT1: STT for object 1 at the final propagation time
% numericSTT2: STT for object 2 at the final propagation time
%
% VARIABLES = 
% mu: Earth gravitational parameter
% rE: Radius of the Earth
% E1: Moment of order 1 for object 1
% E2: Propagated normalised covariance for object 1
% E2_inter: De-normalised final covariance for object 1
% Eq_to_Del1: Equinoctial to Delaunay jacobian for rotation for object 1
% E2_2: E2_inter rotated to Equinoctial frame for object 1
% PQW_to_IJK: PQW to IJK transformation matrix for object 1
% J1: 3x3 Equinoctial to PQW jacobian transformation matrix for object 1
% J2: 6x6 Equinoctial to PQW jacobian transformation matrix for object 1
% E12: Moment of order 1 for object 2
% E22: Propagated normalised covariance for object 2
% E22_inter: De-normalised final covariance for object 2
% Eq_to_Del2: Equinoctial to Delaunay jacobian for rotation for object 2
% E22_2: E2_inter rotated to Equinoctial frame for object 2
% PQW_to_IJK2: PQW to IJK transformation matrix for object 2
% J12: 3x3 Equinoctial to PQW jacobian transformation matrix for object 2
% J22: 6x6 Equinoctial to PQW jacobian transformation matrix for object 2
% 
% OUTPUTS = 
% Q_new: Final propagated covariance for object 1
% Q_new_2: Final propagated covariance for object 2
% =======================================================================

mu = constants.mu;
rE = constants.rE;

E1 = E1calc(constants.STTOrder, Q_bar, mu_old1, constants, numericSTT1);
E2 = E2calc(constants.STTOrder, Q_bar, mu_old1, E1, constants, numericSTT1);
E2_inter = normalize(E2,rE,'mat','Del',0); 
Eq_to_Del1 = computeEqDelPartials(constants, New_plot_COE, 0);
E2_2 = Eq_to_Del1 * E2_inter(1:6,1:6) * Eq_to_Del1';
PQW_to_IJK = PQW_to_IJK_transform(New_plot_COE,1);
J1 = JacobianCalc(New_plot_COE, mu, 0, 1);
J2 = [PQW_to_IJK zeros(3,3); zeros(3,3) PQW_to_IJK];
Q_new = J2 * J1 * E2_2 * J1' * J2';

E12 = E1calc(constants.STTOrder, Q_bar2, mu_old2, constants, numericSTT2);
E22 = E2calc(constants.STTOrder, Q_bar2, mu_old2, E12, constants, numericSTT2);
E22_inter = normalize(E22,rE,'mat','Del',0);
Eq_to_Del2 = computeEqDelPartials(constants, New_plot_COE_2, 0);
E22_2 = Eq_to_Del2 * E22_inter(1:6,1:6) * Eq_to_Del2';
PQW_to_IJK2 = PQW_to_IJK_transform(New_plot_COE_2,1);
J12 = JacobianCalc(New_plot_COE_2, mu, 0, 1);
J22 = [PQW_to_IJK2 zeros(3,3); zeros(3,3) PQW_to_IJK2];
Q_new_2 = J22 * J12 * E22_2 * J12' * J22';

end
