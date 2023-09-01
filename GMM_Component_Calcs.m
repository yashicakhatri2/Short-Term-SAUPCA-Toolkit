function [mean_new, Q_bar, w_bar] = GMM_Component_Calcs(jmax, nu_OG, Q, fig4)
% This function splits the initial Gaussian distribution into smaller GMM
% components.
%
% ================================================================
% INPUTS = 
% jmax: Number of GMM components to split into
% nu_OG: The nominal "mean" or the nominal trajectory
% Q: Nominal epoch covariance
% fig4: Flag to plot weights vs means
% 
% VARIABLES = 
% ao: Semi-major axis from nominal trajectory
% m: GMM component mean distribution constant that spreads means out from Horwood suboptimal algorithm
% sigma: GMM component deviation constant that spreads means out from Horwood suboptimal algorithm
% mu: Initial GMM split means
% M: Optimisation metric 1 from Horwood suboptimal algorithm
% n: Optimisation metric 2 from Horwood suboptimal algorithm
% lb: Lower bounds of the weights
% Aeq: Optimisation parameter = Aeq(1) x1 + Aeq(2) x2 + Aeq(3) x3 + ... = beq
% beq: Optimisation parameter = Aeq(1) x1 + Aeq(2) x2 + Aeq(3) x3 + ... = beq
% optimOptions: Optimisation tolerances 
% w: Optimised weights
% w_tilde: Intermediate weights to acheieve final weights for GMM components - from Horwood suboptimal algorithm
% mu_tilde: Intermediate means to acheieve final means for GMM components - from Horwood suboptimal algorithm
% mu_cap: Intermediate means #2 to acheieve final means for GMM components - from Horwood suboptimal algorithm
% sigma_tilde_2: Intermediate sigma to  acheieve final parameters for GMM components - from Horwood suboptimal algorithm
% w_cap: Intermediate weights #2 to acheieve final weights for GMM components - from Horwood suboptimal algorithm
% sigma_cap_2: Intermediate sigma #2 to  acheieve final parameters for GMM components - from Horwood suboptimal algorithm
% e1: Unit vector in direction 1 
% nu_bar: Intermediate means #3 to acheieve final means for GMM components - from Horwood suboptimal algorithm
% mean_bar: Intermediate means #4 to acheieve final means for GMM components - from Horwood suboptimal algorithm
%
% OUTPUTS = 
% mean_new: Final GMM component means 
% Q_bar: Final GMM component covariance 
% w_bar: Final GMM component weights
% ================================================================

ao = nu_OG(1);
m = 6;
sigma = 2 * m / (jmax - 1);

if sigma < 0 || sigma > 1
    fprintf('Invalid sigma value.\n')
    return;
end

mu = NaN(jmax,1);
M = NaN(jmax,jmax);
n = NaN(jmax,1);
for alpha = 1:jmax
    mu(alpha) = -m + sigma * (alpha - 1);
end
for alpha = 1:jmax
    for beta = 1:jmax
        M(alpha,beta) = normpdf(mu(alpha) - mu(beta),0,sqrt(2 * sigma^2));
    end
    n(alpha) = normpdf(mu(alpha),0,sqrt(sigma^2+1));
end

lb = zeros(jmax,1);
Aeq = ones(1,jmax);
beq = 1;
optimOptions = optimoptions(@quadprog,'ConstraintTolerance', 1E-25,'OptimalityTolerance',1E-25);
w = quadprog(M,-n,[],[],Aeq,beq,lb,[],[], optimOptions);

for alpha = 1:jmax
    w_tilde(alpha) = sqrt(2 * pi) / sqrt(1 - sigma^2) * w(alpha) * exp(mu(alpha)^2 / 2 / (1 - sigma^2));
    mu_tilde = mu(alpha) / (1 - sigma^2);
    mu_cap(alpha) = ao + sqrt(Q(1,1)) * mu_tilde;
end
sigma_tilde_2 = sigma^2 / (1 - sigma^2);
w_cap = w_tilde;
sigma_cap_2 = sigma_tilde_2 * Q(1,1);
e1 = [1; zeros(5,1)];
Q_bar = (sigma_cap_2^-1 * (e1 * e1') + (Q\eye(6))) \ eye(6);
nu_bar = NaN(jmax,6);
w_bar = NaN(jmax,1);
mean_bar = NaN(6,jmax);
mean_new = NaN(6,jmax);
for alpha = 1:jmax
    nu_bar(alpha,:) = Q_bar * (sigma_cap_2^-1 * mu_cap(alpha) * e1 + (Q\eye(6)) * nu_OG);
    w_bar(alpha,:) = w_cap(alpha) * normpdf(mu_cap(alpha) - e1' * nu_OG, 0, sqrt(sigma_cap_2 + e1'*Q*e1));
    mean_bar(:,alpha) = nu_bar(alpha,:)' - nu_OG;
    mean_new(:,alpha) = nu_bar(alpha,:)';
end

w_bar = w_bar / sum(w_bar);

% Plot to see the spread of weights and means
if fig4 == 1
    any(w < 0)
    figure(5)
    semilogy(mu, w, 'bs')
    hold on;
    semilogy(mu, w, 'b--')
    grid on
    title(sprintf('Weight and mean distribution for N = %i',jmax))
    xlabel('x')
    ylabel('Weights')
end
end