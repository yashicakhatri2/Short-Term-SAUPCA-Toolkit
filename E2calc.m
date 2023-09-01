function E2 = E2calc(M, P, m, Edm, constants, numericSTT)
% Compute propagated covariances using STTs.
% =======================================================================
% INPUTS = 
% M: STT order
% P: Covariance
% m: Mean
% Edm: Moment of order 1
% constants:
%   StateDim: Dimension of the state = 8
% numericSTT: STT used to map the state to the final time
% 
% VARIABLES:
% sum: Final sum to compute the moment
% sumo: Sum 0 STT mapped moment
% sum1: Sum 1 STT mapped moment
% sum2: Sum 2 STT mapped moment
% sum3: Sum 3 STT mapped moment
% sum4: Sum 4 STT mapped moment
% sum5: Sum 5 STT mapped moment
% sum6: Sum 6 STT mapped moment
% sum7: Sum 7 STT mapped moment
% sum8: Sum 8 STT mapped moment
% sum9: Sum 9 STT mapped moment
% mterm: Mapping calc
% 
% OUTPUT =
% E1: Mapped final covariance
% =======================================================================

E2 = [];

for i = 1:constants.StateDim
    
    for j = 1:constants.StateDim
        % Initialize
        sum = 0;
        if M == 1
                % p = 1, q = 1
                sum1 = 0;
                for k1 = 1:constants.StateDim
                    sumo = 0;
                    for l1 = 1:constants.StateDim
                        mterm = phi_calc(1, i, k1, numericSTT) * phi_calc(1, j, l1, numericSTT) * Initial_Ecalc([m(k1), m(l1)], P(k1,l1), 2);
                        sumo = sumo + mterm;
                    end
                    sum1 = sum1 + sumo;
                end
                sum = sum + sum1 - Edm(i) * Edm(j);
    % M == 2
        elseif M == 2
                % p = 1, q = 1
                sum1 = 0;
                for k1 = 1:constants.StateDim
                    sumo = 0;
                    for l1 = 1:constants.StateDim
                        mterm = phi_calc(1, i, k1, numericSTT) * phi_calc(1, j, l1, numericSTT) *...
                            Initial_Ecalc([m(k1), m(l1)], P(k1,l1), 2);
                        sumo = sumo + mterm;
                    end
                    sum1 = sum1 + sumo;
                end
                sum = sum + sum1;
                
                % p = 1, q = 2
                sum2 = 0;
                for k1 = 1:constants.StateDim
                    sum1 = 0;
                    for l1 = 1:constants.StateDim
                        sumo = 0;
                        for l2 = 1:constants.StateDim
                            mterm = phi_calc(1, i, k1, numericSTT) * phi_calc(2, j, [l1;l2], numericSTT) *...
                                Initial_Ecalc([m(k1), m(l1), m(l2)], [P(l1,l2), P(k1,l2), P(k1,l1)], 3);
                            sumo = sumo + mterm;
                        end
                        sum1 = sum1 + sumo;
                    end
                    sum2 = sum2 + sum1;
                end
                sum = sum + sum2/2;
                
                % p = 2, q = 1
                sum3 = 0;
                for k1 = 1:constants.StateDim
                    sum1 = 0;
                    for k2 = 1:constants.StateDim
                        sumo = 0;
                        for l1 = 1:constants.StateDim
                            mterm = phi_calc(1, j, l1, numericSTT) * phi_calc(2, i, [k1;k2], numericSTT) *...
                                Initial_Ecalc([m(k1), m(k2), m(l1)], [P(k2,l1), P(k1,l1), P(k1,k2)], 3);
                            sumo = sumo + mterm;
                        end
                        sum1 = sum1 + sumo;
                    end
                    sum3 = sum3 + sum1;
                end
                sum = sum + sum3/2;
                
                % p = 2, q = 2
                sum4 = 0;
                for k1 = 1:constants.StateDim
                    sum1 = 0;
                    for k2 = 1:constants.StateDim
                        sumo = 0;
                        for l1 = 1:constants.StateDim
                            sumo1 = 0;
                            for l2 = 1:constants.StateDim
                                mterm = phi_calc(2, i, [k1;k2], numericSTT) * phi_calc(2, j, [l1;l2], numericSTT) *...
                                    Initial_Ecalc([m(k1), m(k2), m(l1), m(l2)], [P(l1,l2), P(k2,l2), P(k1,l2), P(k2,l1), P(k1,l1), P(k1,k2)], 4);
                                sumo1 = sumo1 + mterm;
                            end
                            sumo = sumo + sumo1;
                        end
                        sum1 = sum1 + sumo;
                    end
                    sum4 = sum4 + sum1;
                end
                sum = sum + sum4/4 - Edm(i) * Edm(j);
    
    
        % M == 3
        elseif M == 3
            % p = 1, q = 1
            sum1 = 0;
            for k1 = 1:constants.StateDim
                sumo = 0;
                for l1 = 1:constants.StateDim
                    mterm = phi_calc(1, i, k1, numericSTT) * phi_calc(1, j, l1, numericSTT) *...
                        Initial_Ecalc([m(k1), m(l1)], P(k1,l1), 2);
                    sumo = sumo + mterm;
                end
                sum1 = sum1 + sumo;
            end
            sum = sum + sum1;
            
            % p = 1, q = 2
            sum2 = 0;
            for k1 = 1:constants.StateDim
                sum1 = 0;
                for l1 = 1:constants.StateDim
                    sumo = 0;
                    for l2 = 1:constants.StateDim
                        mterm = phi_calc(1, i, k1, numericSTT) * phi_calc(2, j, [l1;l2], numericSTT) *...
                            Initial_Ecalc([m(k1), m(l1), m(l2)], [P(l1,l2), P(k1,l2), P(k1,l1)], 3);
                        sumo = sumo + mterm;
                    end
                    sum1 = sum1 + sumo;
                end
                sum2 = sum2 + sum1;
            end
            sum = sum + sum2/2;
            
            % p = 2, q = 1
            sum3 = 0;
            for k1 = 1:constants.StateDim
                sum1 = 0;
                for k2 = 1:constants.StateDim
                    sumo = 0;
                    for l1 = 1:constants.StateDim
                         mterm = phi_calc(1, j, l1, numericSTT) * phi_calc(2, i, [k1;k2], numericSTT) *...
                            Initial_Ecalc([m(k1), m(k2), m(l1)], [P(k2,l1), P(k1,l1), P(k1,k2)], 3);
                        sumo = sumo + mterm;
                    end
                    sum1 = sum1 + sumo;
                end
                sum3 = sum3 + sum1;
            end
            sum = sum + sum3/2;
            
            % p = 2, q = 2
            sum4 = 0;
            for k1 = 1:constants.StateDim
                sum1 = 0;
                for k2 = 1:constants.StateDim
                    sumo = 0;
                    for l1 = 1:constants.StateDim
                        sumo1 = 0;
                        for l2 = 1:constants.StateDim
                            mterm = phi_calc(2, i, [k1;k2], numericSTT) * phi_calc(2, j, [l1;l2], numericSTT) *...
                                Initial_Ecalc([m(k1), m(k2), m(l1), m(l2)], [P(l1,l2), P(k2,l2), P(k1,l2), P(k2,l1), P(k1,l1), P(k1,k2)], 4);
                            sumo1 = sumo1 + mterm;
                        end
                        sumo = sumo + sumo1;
                    end
                    sum1 = sum1 + sumo;
                end
                sum4 = sum4 + sum1;
            end
            sum = sum + sum4/4;

            % p = 1, q = 3
            sum5 = 0;
            for k1 = 1:constants.StateDim
                sum1 = 0;
                for l1 = 1:constants.StateDim
                    sumo = 0;
                    for l2 = 1:constants.StateDim
                        sumo1 = 0;
                        for l3 = 1:constants.StateDim
                            mterm = phi_calc(1, i, k1, numericSTT) * phi_calc(3, j, [l1;l2;l3], numericSTT) *...
                                Initial_Ecalc([m(k1), m(l1), m(l2), m(l3)], [P(l2,l3), P(l1,l3), P(k1,l3), P(l1,l2), P(k1,l2), P(k1,l1)], 4);
                            sumo1 = sumo1 + mterm;
                        end
                        sumo = sumo + sumo1;
                    end
                    sum1 = sum1 + sumo;
                end
                sum5 = sum5 + sum1;
            end
            sum = sum + sum5 / factorial(3);
            
            % p = 3, q = 1
            sum6 = 0;
            for k1 = 1:constants.StateDim
                sum1 = 0;
                for k2 = 1:constants.StateDim
                    sumo1 = 0;
                    for k3 = 1
                        sumo = 0;
                        for l1 = 1:constants.StateDim
                            mterm = phi_calc(1, j, l1, numericSTT) * phi_calc(3, i, [k1;k2;k3], numericSTT) *...
                                Initial_Ecalc([m(k1), m(k2), m(k3), m(l1)], [P(k3,l1), P(k2,l1), P(k1,l1), P(k2,k3), P(k1,k3), P(k1,k2)], 4);
                            sumo = sumo + mterm;
                        end
                        sumo1 = sumo1 + sumo;
                    end
                    sum1 = sum1 + sumo1;
                end
                sum6 = sum6 + sum1;
            end
            sum = sum + sum6 / factorial(3);
            
            % p = 2, q = 3
            sum7 = 0;
            for k1 = 1:constants.StateDim
                sum1 = 0;
                for k2 = 1:constants.StateDim
                    sumo = 0;
                    for l1 = 1:constants.StateDim
                        sumo1 = 0;
                        for l2 = 1:constants.StateDim
                            sumo2 = 0;
                            for l3 = 1:constants.StateDim
                                mterm = phi_calc(2, i, [k1;k2], numericSTT) * phi_calc(3, j, [l1;l2;l3], numericSTT) *...
                                    Initial_Ecalc([m(k1), m(k2) m(l1), m(l2), m(l3)], [P(l2,l3), P(l1,l3), P(l1,l2), P(k2, l3), P(k2,l2), P(k2,l1), P(k1,l3), P(k1,l2), P(k1,l1), P(k1,k2)], 5);
                                sumo2 = sumo2 + mterm;
                            end
                            sumo1 = sumo1 + sumo2;
                        end
                        sumo = sumo + sumo1;
                    end
                    sum1 = sum1 + sumo;
                end
                sum7 = sum7 + sum1;
            end
            sum = sum + sum7 / factorial(2) / factorial(3);
            
            % p = 3, q = 2
            sum8 = 0;
            for k1 = 1:constants.StateDim
                sum1 = 0;
                for k2 = 1:constants.StateDim
                    sumo = 0;
                    for k3 = 1:constants.StateDim
                        sumo1 = 0;
                        for l1 = 1:constants.StateDim
                            sumo2 = 0;
                            for l2 = 1:constants.StateDim
                                mterm = phi_calc(3, i, [k1;k2; k3], numericSTT) * phi_calc(2, j, [l1;l2], numericSTT) *...
                                    Initial_Ecalc([m(k1), m(k2) m(k3), m(l1), m(l2)], [P(l1,l2), P(k3,l2), P(k3,l1), P(k2,l2), P(k2,l1), P(k2,k3), P(k1,l2), P(k1,l1), P(k1,k3), P(k1,k2)], 5);
                                sumo2 = sumo2 + mterm;
                            end
                            sumo1 = sumo1 + sumo2;
                        end
                        sumo = sumo + sumo1;
                    end
                    sum1 = sum1 + sumo;
                end
                sum8 = sum8 + sum1;
            end
            sum = sum + sum8 / factorial(3) / factorial(2);
            
            % p = 3, q = 3
            sum9 = 0;
            for k1 = 1:constants.StateDim
                sum1 = 0;
                for k2 = 1:constants.StateDim
                    sumo = 0;
                    for k3 = 1:constants.StateDim
                        sumo1 = 0;
                        for l1 = 1:constants.StateDim
                            sumo2 = 0;
                            for l2 = 1:constants.StateDim
                                sumo3 = 0;
                                for l3 = 1:constants.StateDim
                                    mterm = phi_calc(3, i, [k1;k2; k3], numericSTT) * phi_calc(3, j, [l1;l2;l3], numericSTT) *...
                                        Initial_Ecalc([m(k1), m(k2) m(k3), m(l1), m(l2), m(l3)], [P(l2,l3), P(l1,l3), P(l1,l2), P(k3,l3), P(k3,l2), P(k3,l1), P(k2,l3), P(k2,l2), P(k2,l1), P(k2,k3), P(k1,l3), P(k1,l2), P(k1,l1), P(k1,k3), P(k1,k2)], 6);
                                    sumo3 = sumo3 + mterm;
                                end
                                sumo2 = sumo2 + sumo3;
                            end
                            sumo1 = sumo1 + sumo2;
                        end
                        sumo = sumo + sumo1;
                    end
                    sum1 = sum1 + sumo;
                end
                sum9 = sum9 + sum1;
            end
            sum = sum + sum9 / factorial(3) / factorial(3) - Edm(i) * Edm(j);
        else
            sum = inf;
            print("STT order exceeds capability")
        end
        E2(i,j) = sum;
        
    end
end
end