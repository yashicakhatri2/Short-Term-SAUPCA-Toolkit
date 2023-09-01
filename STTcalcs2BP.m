function dx_new = STTcalcs2BP(dx, p, constants, numericSTT)
% Mean propagation with STTs.
% =======================================================================
% INPUTS = 
% dx: State to be mapped
% p: STT order
% constants:
%   StateDim: Dimension of the state = 8
% numericSTT: STT used to map the state to the final time
% 
% VARIABLES:
% sum: Final sum to compute the state
% sum1: Sum 1 STT mapped state
% xterm: phi * state = inidividual order mapped state
% sum2: Sum 2 STT mapped state
% 
% OUTPUT =
% dx_new: Mapped final state
% =======================================================================

dx_new = [];

for i = 1:constants.StateDim
    % Initialize
    sum = 0;
    if p == 1
        % 1
        sum1 = 0;
        for k1 = 1:constants.StateDim
            xterm = phi_calc(1, i, k1, numericSTT) * dx(k1);
            sum1 = sum1 + xterm;
        end
        sum = sum + sum1;
        
    elseif p == 2
        % 1
        sum1 = 0;
        for k1 = 1:constants.StateDim
            xterm = phi_calc(1, i, k1, numericSTT) * dx(k1);
            sum1 = sum1 + xterm;
        end
        sum = sum + sum1;

        % 2
        sum2 = 0;
        for k1 = 1:constants.StateDim
            sum1 = 0;
            for k2 = 1:constants.StateDim
               xterm = phi_calc(2, i, [k1;k2], numericSTT) * dx(k1) * dx(k2);
               sum1 = sum1 + xterm;
            end
            sum2 = sum2 + sum1;
        end
        sum = sum + sum2/2;
    end
    dx_new(i) = sum;
end
end