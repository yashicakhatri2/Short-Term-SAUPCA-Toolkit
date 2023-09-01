function OUT = COE_to_Cartesian(IN,mu,dir)
% Convert Classical Orbital Elements to a cartesian state and vice-versa.
% =======================================================================
% INPUTS = 
% IN: Input set - COE or cartesian state
% mu: Gravitational parameter
% dir: 1 - COE to cartesian state, 0 - Cartesian state to COE
% 
% VARIABLES = 
% a: Semi-major axis 
% e: Eccentricity
% i: Inclination
% O: Right Ascension of the Ascending Node
% w: Arguemnt of periapsis
% M: Mean Anomaly
% TA: True Anomaly
% p: Semi-latus rectum
% rPQW: Position vector in PQW transition frame
% vPQW: Velocity vector in PQW transition frame
% PQWtoIJK: PQW to IJK rotation matrix
% rIJK: Position vector in IJK frame
% vIJK: Velocity vector in IJK frame
% State: Cartesian state
% rvec: Position vector from input
% r: Position norm from input
% vvec: Velocity vector from input
% v: Velocity norm from input
% hvec: Angular momentum vector
% h: Angular momentum vector norm
% K: Vector in khat direction [0 0 1]
% nvec: Cross vector between K and angular momentum vector
% evec: Eccentricity vector
% E: Eccentric anomaly
% 
% OUTPUT =
% OUT: Output set - COE or cartesian state
% =======================================================================

if dir == 1    
    a = IN(1); e = IN(2); i = IN(3); O = IN(4); w = IN(5); M = IN(6);
    [~,TA] = MAtoEccTA(M,e);
    p = a * (1 - e^2);
    if abs(imag(w)) > 0 || abs(imag(O)) > 0 || abs(imag(M)) > 0
        return
    end
    
    rPQW = [p * cos(TA)/(1+e*cos(TA)); p * sin(TA)/(1+e*cos(TA)); 0];
    vPQW = [-sqrt(mu/p)*sin(TA); sqrt(mu/p)*(e+cos(TA)); 0];
    PQWtoIJK = PQW_to_IJK_transform(IN, 1);
    rIJK = PQWtoIJK * rPQW;
    vIJK = PQWtoIJK * vPQW;

    OUT = [rIJK; vIJK];
else
    State = IN;
    rvec = State(1:3); r = norm(rvec);
    vvec = State(4:6); v = norm(vvec);
    hvec = cross(rvec, vvec); h = norm(hvec); K = [0;0;1];
    nvec = cross(K,hvec);
    evec = ((v^2 - mu/r) * rvec - (dot(rvec,vvec)) * vvec) / mu; e = norm(evec);
    E = v^2/2 - mu/r;
    if e ~= 1
        a = -mu/2/E;
        p = a * (1 - e^2);
    else
        p = h^2/mu;
        a = inf;
    end
    i = acosd(hvec(3)/h);
    O = acosd(nvec(1)/norm(nvec));
    if nvec(2) < 0
        O = 360 - O;
    end
    if abs(dot(nvec, evec)/norm(nvec)/e) > 0.9999999999999 && abs(dot(nvec, evec)/norm(nvec)/e) < 1.00000001
        w = 0;
    else
        w = acosd(dot(nvec, evec)/norm(nvec)/e);
        if evec(3) < 0
            w = 360 - w;
        end
    end
    if abs(dot(evec,rvec)/e/r) > 1
        TA = acosd(1 * sign(dot(evec,rvec)/e/r));
    else
        TA = acosd(dot(evec,rvec)/e/r);
    end
    if dot(rvec,vvec) < 0
        TA = 360 - TA;
    end
    E = 2 * atan(sqrt((1-e)/(1+e)) * tand(TA/2));
    M = E - e * sin(E);
    if abs(imag(w)) > 0 || abs(imag(O)) > 0 || abs(imag(M)) > 0
        return
    end
    OUT = [a; e; i * pi / 180; wrapToPi(O * pi / 180); wrapToPi(w * pi / 180); wrapToPi(M)];
end
end