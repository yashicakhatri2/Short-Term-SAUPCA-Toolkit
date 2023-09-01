function Cart = Equinoctial_to_Cartesian(Eq, mu)
% Convert Equinoctial set to a cartesian state.
% =======================================================================
% INPUTS = 
% Eq: Equinoctial element set
% mu: Gravitational parameter
% 
% SUB-FUNCTIONS = 
% COE_to_Equinoctial: Converts Equinoctial element set to Classical Orbital Element set
% COE_to_Cartesian: Converts Classical Orbital Element set to Cartesian state
%
% OUTPUT =
% Cart: Output Cartesian state
% =======================================================================

    COE =  COE_to_Equinoctial(Eq, mu, 0);
    Cart = COE_to_Cartesian(COE, mu, 1);
end