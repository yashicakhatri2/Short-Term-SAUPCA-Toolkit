function [COE2, P_prop] = define_cases(case_flag)
% Defined test cases that have pre-saved results to confirm functionality.
% =======================================================================
% INPUTS = 
% case_flag: Test number from define_cases function 
% 
% OUTPUTS = 
% COE2: Classical orbital element set for object 2 in conjunction
% P_prop: Nominal Time of Closest Approach for the two object in conjunction
% =======================================================================

if case_flag == 1
    COE2 = [9843.33683736869         0.241038139941096          1.07980901556839        -0.191768531625329         0.970045507875877         -1.52824344879098]';
    P_prop = 1.5 * 24 * 3600;
end

end