function [gradLogGaussVA] = grad_GaussVA(VAMu, VAB, VAD, parameter)
% Compute gradients of Gaussian Factor Variational Distribution
    [P, K] = size(VAB);
    
    inv_D2 = diag(VAD.^(-2));

    A = parameter - VAMu;
    gradLogGaussVA = - ( inv_D2*A - inv_D2*VAB*( eye(K) + VAB'* inv_D2 *VAB )^(-1)*VAB'*inv_D2 * A );
end