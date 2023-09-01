function PQW_to_IJK = PQW_to_IJK_transform(COE, dir)
% This function calculates the rotation matrix between PQW and IJK frames.
%
% =======================================================================
% INPUTS = 
% COE: Current classical orbital element set
% dir: 1 - PQW to IJK rotation, 0 - IJK to PQW rotation
% 
% VARIABLES = 
% i: Inclination
% O: Right Ascension of the Ascending Node
% w: Arguemnt of periapsis
% R1: Rotation 1
% R2: Rotation 2
% R3: Rotation 3
% 
% OUTPUTS = 
% PQW_to_IJK: Rotation matrix between PQW and IJK, depending on the direction input
% =======================================================================

i = COE(3); O = COE(4); w = COE(5);

R1 = [cos(-O) sin(-O) 0;...
    -sin(-O) cos(-O) 0;
    0 0 1];
R2 = [1 0 0;...
    0 cos(-i) sin(-i);...
    0 -sin(-i) cos(-i)];
R3 = [cos(-w) sin(-w) 0;...
    -sin(-w) cos(-w) 0;...
    0 0 1];

if dir == 1
    PQW_to_IJK = R1 * R2 * R3;
else
    PQW_to_IJK = (R1 * R2 * R3)\eye(3);
end


end