function Out = COE_to_Equinoctial(COE, mu, flag)
% This function converts COE to Equinoctial state, and vice-versa.
%
% =======================================================================
% INPUTS = 
% COE: Current classical orbital element set
% mu: Gravitational parameter of the Earth
% flag: 1 - Convert from COE to Equinoctial state, 0 - Convert from Equinoctial to COE state
% 
% VARIABLES = 
% a: Semi-major axis  
% e: Eccentricity
% i: Inclination
% O: Right Ascension of the Ascending Node
% w: Arguemnt of periapsis
% M: Mean Anomaly at time, t
% h: Equinoctial element 2
% k: Equinoctial element 3
% p: Equinoctial element 4
% q: Equinoctial element 5
% lam: Equinoctial element 6
% 
% OUTPUTS = 
% Out: Output state in COE or Equinoctial state, depending on flag
% =======================================================================

if flag == 1
    a = COE(1); e = COE(2); i = COE(3); O = COE(4); w = COE(5); M = COE(6); 
    h = (e * sin(w + O));
    k = (e * cos(w + O));
    p = (tan(i/2) * sin(O));
    q = (tan(i/2) * cos(O));
    if abs(imag(w)) > 0 || abs(imag(O)) > 0 || abs(imag(M)) > 0 || e < 0
        return
    end
    lam = wrapToPi(M + w + O);
    Out = [a; h; k; p; q; lam];
else
    a = COE(1); h = COE(2); k = COE(3); p = COE(4); q = COE(5); lam = COE(6);
    O = wrapToPi(atan2(p,q));
    w = wrapToPi(atan2(h,k) - O);
    if w + O == 0
        e = k / cos(w + O);
    else
        e = h / sin(w + O);
    end
    
    if p == 0
       i = wrapToPi(2 * atan2(q,cos(O)));
    else
       i = wrapToPi(2 * atan2(p,sin(O)));
    end
    
    M = wrapToPi(lam - w - O);
    
    if abs(imag(lam)) > 0 || abs(imag(w)) > 0 || abs(imag(O)) > 0 || abs(imag(M)) > 0 || e < 0
        return
    end
    
    Out = [a; e; i; O; w; M];
end

end