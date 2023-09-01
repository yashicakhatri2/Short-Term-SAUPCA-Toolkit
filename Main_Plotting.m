function Main_Plotting(results,constants,GMM_mean,Q_new)
% This function calculates the Jacobian between the Equinoctial element set
% and the PQW frame.
%
% =======================================================================
% INPUTS = 
% constants: Initial defined constants
%   points: Number of points initially defined to compute MC_Pc  
%   case_flag: Test number from define_cases function, using predefined test cases
%   JMAX: List of predefined number of GMMs components for each object used to compute GMMSTT_Pc
%   testFig2: Plot - Testing the GMM spread at a set time t = P_prop from function define cases. JMAX FIXED AT 15 COMPONENTS. DOES NOT PROVIDE ACCURATE PROBABILITY OF COLLISION RESULTS.
%   plot_GMM_ell: 1 - plot GMM ellipses (usually used when testFlagFig2 = 1 and testing at a set time t)
% results:
%   MC_Pc: Monte Carlo probability of collision
%   points: Number of points used to compute MC_Pc
%   MCRunTime: Runtime for MC_Pc calculations, some are saved, use tic-toc to compute and save new ones
%   GMMSTT_Pc: Semi-analytical probability of collision
%   JMAX: List of number of GMMs components for each object used to compute GMMSTT_Pc
%   GMMSTTRunTime: List of times to compute GMMSTT_Pc using # components listed in JMAX
% GMM_MEAN: GMM component propagated means  
% Q_new: GMM component propagated covariance
% 
% VARIABLES = 
% MC_Pc: Monte Carlo probability of collision
% points: Number of points used to compute MC_Pc
% JMAX: List of number of GMMs components for each object used to compute GMMSTT_Pc
% PC_GMMi: Semi-analytical probability of collision
% case_flag: Test number from define_cases function, using predefined test cases
% fig2: Plot - Testing the GMM spread at a set time t = P_prop from function define cases. JMAX FIXED AT 15 COMPONENTS. DOES NOT PROVIDE ACCURATE PROBABILITY OF COLLISION RESULTS.
% plot_GMM_ell: 1 - plot GMM ellipses (usually used when testFlagFig2 = 1 and testing at a set time t)
% jmax: Last component of list of number of GMMs components for each object used to compute GMMSTT_Pc
% MSPColor: Color 1 for plotting
% CartColor: Color 2 for plotting
% Other plotting variables
% 
% OUTPUTS = 
% J: Desired plots
% =======================================================================


MC_Pc = results.MC_Pc;
points = results.points;
JMAX = results.JMAX;
PC_GMMi = results.GMMSTT_Pc;
case_flag = constants.case_flag;
fig2 = constants.testFig2;
plot_GMM_ell = constants.plot_GMM_ell;
jmax = results.JMAX(end);

MSPColor = '#EF1075';
CartColor = '#1075ef';

if fig2 ~= 1
    MC_Std = 1.960 * sqrt(MC_Pc * (1 - MC_Pc) / points);
    
    gcf = figure(1);
    str = sprintf('%.2fE%i SDS+SP Monte Carlo Points Pc',points / 10^floor(log10(points)),floor(log10(points)));
    CI1 = MC_Pc + MC_Std;
    CI2 = MC_Pc - MC_Std;
    semilogy(0:JMAX(end) + 20, ones(1,JMAX(end)+21) * MC_Pc,'Color',CartColor, 'LineWidth', 2)
    hold on
    plot(0:JMAX(end) + 20, ones(1,JMAX(end)+21) * CI1, ':','Color',CartColor,  'LineWidth', 1.2)
    plot(JMAX,PC_GMMi,'d','LineWidth',1,'MarkerSize',7,'markerfacecolor',MSPColor,'MarkerEdgeColor','k')
    plot(0:JMAX(end) + 20, ones(1,JMAX(end)+21) * CI2, ':','Color',CartColor,  'LineWidth', 1.2)
    legend(str,'95% CI Bounds','SDS+SP GMM-STT Method Pc')
    xlabel('# GMM Mixtures'); ylabel('Probability of Collision')
    grid on;
    xlim([-inf, inf])
    if results.MC_Pc ~= 0
        ylim([results.MC_Pc*0.5,results.MC_Pc*1.5])
    end
    set(gcf, 'Position',  [200, 150, 700, 250])
    
    PC_GMMi_Rel = abs(PC_GMMi - MC_Pc) ./ MC_Pc * 100;
    MaxRelError = (MC_Std/MC_Pc) * 100;
    
    gcf4 = figure(4);
    plot(0:JMAX(end) + 20, ones(1,JMAX(end)+21) * MaxRelError, '--','Color',CartColor,  'LineWidth', 2)
    hold on
    plot(JMAX,PC_GMMi_Rel,'d','LineWidth',1,'MarkerSize',7,'markerfacecolor',MSPColor,'MarkerEdgeColor','k')
    grid on
    xlim([-inf, inf])
    if MC_Pc ~= 0
        ylim([0,MaxRelError * 1.2])
    end
    xlabel('# GMM Mixtures'); ylabel('Relative Error [%]')
    title('Relative Error Threshold')
    legend('Relative Error for 95% CI Bounds','Current Relative Error')
    set(gcf4, 'Position',  [200, 150, 700, 250])
end
if fig2 == 1
% Ellipsoid calculations
    syms x y
    if plot_GMM_ell == 1
        for i = 1:jmax
            Pxx = Q_new{i}(1,1);
            Pyy = Q_new{i}(2,2);
            Pzz = Q_new{i}(3,3);
            Pxy = Q_new{i}(1,2);
            Pyz = Q_new{i}(2,3);
            Pxz = Q_new{i}(1,3);
            xoff = GMM_mean(1,i);
            yoff = GMM_mean(2,i);
            zoff = GMM_mean(3,i);
            fhxy{i} = @(x,y) (Pyy .* (x-xoff)^2 - 2 .* Pxy .* (x-xoff) .* (y-yoff) + Pxx .* (y-yoff)^2)./(Pxx * Pyy - Pxy^2) - chi2inv(.997,1);
            fhxz{i} = @(x,z) (Pzz .* (x-xoff)^2 - 2 .* Pxz .* (x-xoff) .* (z-zoff) + Pxx .* (z-zoff)^2)./(Pxx * Pzz - Pxz^2) - chi2inv(.997,1);
            fhyz{i} = @(y,z) (Pzz .* (y-yoff)^2 - 2 .* Pyz .* (y-yoff) .* (z-zoff) + Pyy .* (z-zoff)^2)./(Pyy * Pzz - Pyz^2) - chi2inv(.997,1);
        end
    end
    
    if case_flag == 1 && constants.usePredefinedCases == 1
        lim1 = [-2000,-600,3990,4050];
        lim2 = [-2000,-600,6780,6830];
        lim3 = [3990,4050,6780,6830];
    else
        fprintf('*******Need to define ellipse plotting upper and lower limits in Main_Plotting.m to plot ellipses.*********\n')
        plot_GMM_ell = 0;
    end

    figure(3)
    subplot(1,3,1)
    xlabel('x')
    ylabel('y')
    plot(GMM_mean(1,:), GMM_mean(2,:), 'k+', 'MarkerSize', 7, 'LineWidth', 1); hold on;
    if plot_GMM_ell == 1
        fig1 = ezplot(fhxy{i}, lim1);
        set(fig1,'color','red','LineStyle', '-','LineWidth',0.5);
        hold on;
        for i = 2:jmax
            fig1 = ezplot(fhxy{i}, lim1);
            set(fig1,'color','red','LineStyle', '-','LineWidth',0.5);
        end
    end
    grid on;
    title('x-y')
    xlim([-inf,inf])
    ylim([-inf,inf])
    
    subplot(1,3,2)
    plot(GMM_mean(1,:), GMM_mean(3,:), 'k+', 'MarkerSize', 7, 'LineWidth', 1); hold on;
    xlabel('x')
    ylabel('z')
    if plot_GMM_ell == 1
        ylim([-200,200])
        fig1 = ezplot(fhxz{i}, lim2);
        set(fig1,'color','red','LineStyle', '-','LineWidth',0.5);
        hold on;
        for i = 2:jmax
            fig1 = ezplot(fhxz{i}, lim2);
            set(fig1,'color','red','LineStyle', '-','LineWidth',0.5);
        end
    end
    grid on;
    title('x-z')
    xlim([-inf,inf])
    ylim([-inf,inf])

    subplot(1,3,3)
    plot(GMM_mean(2,:), GMM_mean(3,:), 'k+', 'MarkerSize', 7, 'LineWidth', 1); hold on;
    xlabel('y')
    ylabel('z')
    if plot_GMM_ell == 1
        fig1 = ezplot(fhyz{i}, lim3);
        set(fig1,'color','red','LineStyle', '-','LineWidth',0.5);
        hold on;
        for i = 2:jmax
            fig1 = ezplot(fhyz{i}, lim3);
            set(fig1,'color','red','LineStyle', '-','LineWidth',0.5);
        end
    end
    
    grid on;
    title('y-z')
    xlim([-inf,inf])
    ylim([-inf,inf])

end
end