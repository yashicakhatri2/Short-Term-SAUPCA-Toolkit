function E1 = E1calc(p, P, m, constants, numericSTT)
% Compute propagated means using STTs.
% =======================================================================
% INPUTS = 
% p: STT order
% P: Covariance
% m: Mean
% constants:
%   StateDim: Dimension of the state = 8
% numericSTT: STT used to map the state to the final time
% 
% VARIABLES:
% sum: Final sum to compute the state
% sum1: Sum 1 STT mapped state
% sum2: Sum 2 STT mapped state
% sum3: Sum 3 STT mapped state
% mterm: Mapping calc
% 
% OUTPUT =
% E1: Mapped final mean
% =======================================================================

E1 = [];
P_dummy = [];

for i = 1:constants.StateDim
    % Initialize
    sum = 0;
    if p == 1
        % 1
        sum1 = 0;
        for k1 = 1:constants.StateDim
            mterm = phi_calc(1, i, k1, numericSTT) * Initial_Ecalc(m(k1), P_dummy, 1);
            sum1 = sum1 + mterm;
        end
        sum = sum + sum1;
        
    elseif p == 2
        % 1
        sum1 = 0;
        for k1 = 1:constants.StateDim
            mterm = phi_calc(1, i, k1, numericSTT) * Initial_Ecalc(m(k1), P_dummy, 1);
            sum1 = sum1 + mterm;
        end
        sum = sum + sum1;
        
        % 2
        sum2 = 0;
        for k1 = 1:constants.StateDim
            sum1 = 0;
            for k2 = 1:constants.StateDim
               mterm = phi_calc(2, i, [k1;k2], numericSTT) * Initial_Ecalc([m(k1), m(k2)], P(k1,k2), 2);
               sum1 = sum1 + mterm;
            end
            sum2 = sum2 + sum1;
        end
        sum = sum + sum2/2;
        
    elseif p == 3
        % 1
        sum1 = 0;
        for k1 = 1:constants.StateDim
            mterm = phi_calc(1, i, k1, numericSTT) * Initial_Ecalc(m(k1), P_dummy, 1);
            sum1 = sum1 + mterm;
        end
        sum = sum + sum1;
        
        % 2
        sum2 = 0;
        for k1 = 1:constants.StateDim
            sum1 = 0;
            for k2 = 1:constants.StateDim
               mterm = phi_calc(2, i, [k1;k2], numericSTT) * Initial_Ecalc([m(k1), m(k2)], P(k1,k2), 2);
               sum1 = sum1 + mterm;
            end
            sum2 = sum2 + sum1;
        end
        sum = sum + sum2/2;
        
        % 3
        sum3 = 0;
        for k1 = 1:constants.StateDim
            sum2 = 0;
            for k2 = 1:constants.StateDim
               sum1 = 0;
               for k3 = 1:constants.StateDim
                   mterm = phi_calc(3, i, [k1;k2; k3], numericSTT) * Initial_Ecalc([m(k1), m(k2), m(k3)], [P(k2,k3), P(k1,k3), P(k1,k2)], 3);
                   sum1 = sum1 + mterm;
               end
               sum2 = sum2 + sum1;
            end
            sum3 = sum3 + sum2;
        end
        sum = sum + sum3 / factorial(3);
  
    else
        sum = inf;
        print("STT order exceeds capability")
    end
    E1(i) = sum;
end
end