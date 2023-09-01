function Pc = Compute_MC_Pc_Cart_1on1(ICState,constants)
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
% P: Cartesian covariance for object 1
% P2: Cartesian covariance for object 2
% mu: Earth gravitaiton parameter 
% State1: Cartesian state for object 1
% State2: Cartesian state for object 2
% points: Number of randomly generated points for each object for the MC analysis
% R_star: Combined hard body radius (r1 + r2)
% P_prop: Nominal time of closest approach (TCA)
% fig2: Plot - Testing the GMM spread at a set time t = P_prop from function define cases. JMAX FIXED AT 15 COMPONENTS.
% n_collision: Number of trajectories in the Monte Carlo analysis that lead to collisions 
% x_new_t_Cartesian: Final cartesian state of object 1 if plotting trajectories at the nominal TCA
% x_new_t_Cartesian_2: Final cartesian state of object 1 if plotting trajectories at the nominal TCA
% counter: Iteration counter
% MC_Points_Cart: Randomly generated cartesian state for object 1
% MC_Points_Cart_2: Randomly generated cartesian state for object 2
% COE1_new: COE set for MC_Points_Cart
% COE2_new: COE set for MC_Points_Cart_2
% S1: Cartesian state of object 1 at the actual Time of Closest Approach (TCA)
% S2: Cartesian state of object 2 at the actual Time of Closest Approach (TCA)
% xe: Distance between the objects at the actual TCA
% Color1: Plotting custom color 1
% Color2: Plotting custom color 2
% 
% OUTPUT =
% Pc: Probability of collision computed using the Monte Carlo method
% =======================================================================

P = ICState.obj1.PoCart;
P2 = ICState.obj2.PoCart;

mu = constants.mu;
State1 = ICState.obj1.Cart;
State2 = ICState.obj2.Cart;
points = constants.points;
R_star = constants.R_star;
P_prop = constants.P_prop;
fig2 = constants.testFig2;

n_collision = 0;
    
x_new_t_Cartesian = NaN(6,points);
x_new_t_Cartesian_2 = NaN(6,points);
counter = 0;

if constants.par == 0
    for i = 1:points
        MC_Points_Cart = mvnrnd(State1, P); 
        MC_Points_Cart_2 = mvnrnd(State2, P2);

        COE1_new = COE_to_Cartesian(MC_Points_Cart,mu,0);
        COE2_new = COE_to_Cartesian(MC_Points_Cart_2,mu,0);
        
        [S1, S2] = FindTCAMCMethod(COE1_new, COE2_new, P_prop, constants);
        if fig2 == 1
            counter = counter + 1;
            x_new_t_Cartesian(:,counter) = S1';
            x_new_t_Cartesian_2(:,counter) = S2';
        end
        xe = norm(S2(1:3) - S1(1:3));

        if xe < R_star
            n_collision = n_collision + 1;
        end
    end
else
    parpool('local',constants.nodes); 
    parfor i = 1:points
        MC_Points_Cart = mvnrnd(State1, P);
        MC_Points_Cart_2 = mvnrnd(State2, P2);

        COE1_new = COE_to_Cartesian(MC_Points_Cart,mu,0);
        COE2_new = COE_to_Cartesian(MC_Points_Cart_2,mu,0);
       
        [S1, S2] = FindTCAMCMethod(COE1_new, COE2_new, P_prop, constants);
        xe = norm(S2(1:3) - S1(1:3));
        if xe < R_star
            n_collision = n_collision + 1;
        end
    end
    delete(gcp('nocreate')) 
end
Pc = sum(n_collision) / points;

%% Plotting

if fig2 == 1
    Color1 = '#1075ef';
    Color2 = '#EF8A10';

    figure(2)
    plot3(x_new_t_Cartesian(1,:),x_new_t_Cartesian(2,:),x_new_t_Cartesian(3,:), '.','Color',Color1,'MarkerSize',10); hold on;
    plot3(x_new_t_Cartesian_2(1,:),x_new_t_Cartesian_2(2,:),x_new_t_Cartesian_2(3,:), '.','Color',Color2,'MarkerSize',10); hold on;
    legend('Object 1 distribution at nominal TCA','Object 2 distribution at nominal TCA')
    grid on;
    title('Isometric View of 2 object distributions')
    xlim([-inf,inf])
    ylim([-inf,inf])
    zlim([-inf,inf])
    xlabel('x (km)')
    ylabel('y (km)')
    zlabel('z (km)')

    figure(3)
    subplot(1,3,1)
    plot(x_new_t_Cartesian(1,:),x_new_t_Cartesian(2,:), '.','Color',Color1); hold on;
    plot(x_new_t_Cartesian_2(1,:),x_new_t_Cartesian_2(2,:), '.','Color',Color2); hold on;

    grid on;
    title('x-y')
    xlim([-inf,inf])
    ylim([-inf,inf])
    xlabel('x')
    ylabel('y')

    subplot(1,3,2)
    plot(x_new_t_Cartesian(1,:),x_new_t_Cartesian(3,:), '.','Color',Color1); hold on;
    plot(x_new_t_Cartesian_2(1,:),x_new_t_Cartesian_2(3,:), '.','Color',Color2); hold on;
    xlabel('x')
    ylabel('z')
    grid on;
    title('x-z')
    xlim([-inf,inf])
    ylim([-inf,inf])

    subplot(1,3,3)
    plot(x_new_t_Cartesian(2,:),x_new_t_Cartesian(3,:), '.','Color',Color1); hold on;
    plot(x_new_t_Cartesian_2(2,:),x_new_t_Cartesian_2(3,:), '.','Color',Color2); hold on;
    xlabel('y')
    ylabel('z')
    grid on;
    title('y-z')

end

end
