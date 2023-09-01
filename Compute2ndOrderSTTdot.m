function phiDotVec = Compute2ndOrderSTTdot(S, AMat2, dim, AMat1, STM1)
% Compute 2nd order STT Equations of Motion (EOMs or STTdot).
% =======================================================================
% INPUTS = 
% S: First order STTdot in a column vector
% AMat2: Linear Dynamics Tensor (LDT) of order 2
% dim: Linear dimension of matrix
% AMat1: Linear Dynamics Tensor (LDT) of order 1
% STM1: State Transition Matrix (STM), which is order 1
% 
% VARIABLES = 
% test1: Transition step to convert S to a 3D tensor
% STM2: S converted to a 3D tensor
% PNewCol: Transition step to convert phiDot to a column vector
% phiDot: Second order STTdot 
% Pmat: 2D part of phiDot used in conversion to a column vector
% 
% OUTPUT =
% phiDotVec: Outputs the STTdot in a column vector form for propagation
% =======================================================================

    test1 = reshape(S,dim^2,dim);
    STM2 = NaN(8,8,8);
    for i = 1:dim
        STM2(:,:,i) = reshape(test1(:,i),dim,dim)';
    end
    PNewCol = NaN(64,8);
    phiDot = pagemtimes(AMat1,STM2) + pagemtimes(AMat2,pagemtimes(STM1,STM1));
    for n = 1:dim
        Pmat(:,:) = phiDot(:,:,n);
        PNewCol(:,n) = reshape(Pmat',dim^2,1);
    end
    phiDotVec = reshape(PNewCol,dim^3,1);
    
end