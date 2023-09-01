function OUT = propagateWithDynamics(COE, STMoi, constants, propTime, flagPropagateSTT)
% This function propagates an input state and/or STT to the desired time 
% using SDS dynamics.
% ========================================================================
% INPUTS = 
% COE: Classical orbital elements of the object to propagate
% STMoi: Previously propagated state + STT to continue propagation
% constants: 
%   STTDim: One dimension of the n-dimensional State Transition Tensor
%   options: Option to define the ode propagation tolerances
%   mu: Earth gravitational parameter
%   rE: Radius of the Earth
%   STTOrder: Order of the STT
% propTime: 
%   (1): Start propagation time
%   (2): End propagation time
% flagPropagateSTT: Flag to propagate STT (1 if yes, 0 if no)
%  
% VARIABLES = 
% dim: One dimension of the n-dimensional State Transition Tensor
% options: Option to define the ode propagation tolerances
% propTime1: Start propagation time
% propTime2: End propagation time
% DelOG: Delaunay form of the input COE set
% Del: Extending the 6 dimensional DelOG to 8 dimensions to incorporate the augmented Hamiltonian set
% Offset: Offset to achieve mean state at initial time for propagation
% Delnorm: Normalised Delaunay state
% STMo: Reshaped STMoi
% PNewCol: Intermediate variable that allows STMoi reshaping
% Pmat: Intermediate variable #2 that allows STMoi reshaping
% STM2Vec: Intermediate variable that allows 2nd order part of STMoi reshaping
% So: Full initial state to start propagation
% tf: Time output from ode propagation
% Sf: State output from ode propagation
% DelfTemp: De-normalised final Delaunay state 
% interOffset: Offset at final time to reverse to osculating state
% Sfend: Final state after implementing offset
% STTf: Final reshaped n^m State Transition Tensor
% test1: Intermediate variable that allows reshaping of the final tensor
% 
% OUTPUTS = 
% OUT1:
%   {1}: Final Delaunay state
%   {2}: Final STT
%   {3}: Final Offset 
% OUT2:
%   {1}: Final Delaunay state
%   {2}: Final offset
% OUT: OUT1 or OUT2 depending on Flag for propagating STT alongside the state
% ========================================================================

% Constants
dim = constants.STTDim;
options = constants.options;
propTime1 = propTime(1);
propTime2 = propTime(2);

% Initial Del Offset calcs
DelOG = COE_to_Delaunay(COE,constants.mu,1);
Del = [DelOG; 0; 0];
if propTime(1) == 0
    Offset = getOffset([DelOG;0;0],0,constants,0);
    Del = [DelOG; 0; 0] + Offset;
end
Delnorm = normalize(Del,constants.rE,'vec','Del',1);

% Propagator
if flagPropagateSTT == 1
    % Initiate STT
    STMo = reshape(STMoi{1}',64,1);
    PNewCol = NaN(64,8);
    if constants.STTOrder == 2
        for n = 1:dim
            Pmat(:,:) = STMoi{2}(:,:,n);
            PNewCol(:,n) = reshape(Pmat',dim^2,1);
        end
        STM2Vec = reshape(PNewCol,dim^3,1);
        STMo = [STMo;STM2Vec];
    end
    
    % Propagate State and STT
    So = [Delnorm;STMo];
    if propTime1 ~= 0
        So = [So(1:8)'-STMoi{3}, So(9:72)', So(73:end)'];
    end
    [tf, Sf] = ode113(@(t,S) SDSDynamics(t,S,constants), [propTime1; propTime2], So, options);
    DelfTemp = normalize(Sf(end,1:8),constants.rE,'vec','Del',0);
    interOffset = normalize(getOffset(DelfTemp,tf(end),constants,1),constants.rE,'vec','Del',1)';
    Sfend = [Sf(end,1:8) + interOffset,Sf(end,9:end)];
    
    % Reshape STT for output
    if constants.STTOrder == 1
        STTf{1} = reshape(Sfend(9:end),dim,dim)';
    elseif constants.STTOrder == 2
        STTf{1} = reshape(Sfend(9:72),dim,dim)';
        test1 = reshape(Sfend(73:end),dim^2,dim);
        for i = 1:dim
            STTf{2}(:,:,i) = reshape(test1(:,i),dim,dim)';
        end
    end
    Delfnorm = Sfend(1:6);
    
    Delf = normalize(Delfnorm,constants.rE,'vec','Del',0);
    OUT1{1} = [wrapToPi(Delf(1:3)) Delf(4:6)];
    OUT1{2} = STTf;
    OUT1{3} = interOffset;
    OUT = OUT1;

else % Propagate State
    So = Delnorm;
    if propTime1 ~= 0
        So = So-STMoi';
    end
    [tf, Sf] = ode113(@(t,S) SDSDynamicsStateOnly(t,S,constants), [propTime1; propTime2], So, options);
    DelfTemp = normalize(Sf(end,1:8),constants.rE,'vec','Del',0);
    interOffset = normalize(getOffset(DelfTemp,tf(end),constants,1),constants.rE,'vec','Del',1)';
    Sfend = Sf(end,1:8) + interOffset;
        
    Delfnorm = Sfend(1:6);
    Delf = normalize(Delfnorm,constants.rE,'vec','Del',0);
    OUT2{1} = [wrapToPi(Delf(1:3)) Delf(4:6)];
    OUT2{2} = interOffset;
    OUT = OUT2;
end

end
