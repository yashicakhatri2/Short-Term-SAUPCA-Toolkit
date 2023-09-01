function OUT = normalize(IN,rE,type1,type2,dir)
% This function allows us to normalise or de-normalise an Equinoctial or
% Delaunay state.
% ========================================================================
% INPUTS = 
% IN: Input in desired format
% rE: Radius of the Earth to normalise distances
% type1: 'vec' - Vector input, 'mat' - matrix input
% type2: 'Equ' - Equinoctial state, 'Del' - Delaunay state
% dir: 1 - normalise, 0 - de-normalise
% 
% VARIABLES = 
% Q_bar: Normalised matrix
% Q: De-normalised matrix
% vec: Normalised or de-normalised vector
%
% OUTPUTS = 
% OUT: Normalised or de-normalised vector or matrix in Equinoctial or Delaunay frame
% ========================================================================

if type2 == 'Equ'
    if type1 == 'mat'
        if dir == 1
            Q = IN;
            Q_bar = Q;
            Q_bar(1,:) = Q(1,:) / rE;
            Q_bar(:,1) = Q(:,1) / rE;
            Q_bar(1,1) = Q(1,1) / rE^2;
            OUT = Q_bar;
        else
            Q_bar = IN;
            Q = Q_bar;
            Q(1,:) = Q_bar(1,:) * rE;
            Q(:,1) = Q_bar(:,1) * rE;
            Q(1,1) = Q_bar(1,1) * rE^2;
            OUT = Q;
        end
    elseif type1 == 'vec'
        vec = IN;
        if dir == 1
            vec(1) = vec(1) / rE;
        else
            vec(1) = vec(1) * rE;
        end
        OUT = vec;
    end
elseif type2 == 'Del'
    if type1 == 'mat'
        if dir == 1
            Q = IN;
            Q_bar = Q;
            Q_bar(4,:) = Q_bar(4,:) / rE^2 .* 3600;
            Q_bar(:,4) = Q_bar(:,4) / rE^2 .* 3600;
            Q_bar(5,:) = Q_bar(5,:) / rE^2 .* 3600;
            Q_bar(:,5) = Q_bar(:,5) / rE^2 .* 3600;
            Q_bar(6,:) = Q_bar(6,:) / rE^2 .* 3600;
            Q_bar(:,6) = Q_bar(:,6) / rE^2 .* 3600;
            OUT = Q_bar;
        else
            Q_bar = IN;
            Q = Q_bar;
            Q(4,:) = Q(4,:) * rE^2 ./3600;
            Q(:,4) = Q(:,4) * rE^2 ./3600;
            Q(5,:) = Q(5,:) * rE^2 ./3600;
            Q(:,5) = Q(:,5) * rE^2 ./3600;
            Q(6,:) = Q(6,:) * rE^2 ./3600; 
            Q(:,6) = Q(:,6) * rE^2 ./3600;
            OUT = Q;
        end
    elseif type1 == 'vec'
        vec = IN;
        if dir == 1
            vec(4:6) = vec(4:6) ./ rE^2 .* 3600;
        else
            vec(4:6) = vec(4:6) .* rE^2 ./3600;
        end
        OUT = vec;
    end
end
end