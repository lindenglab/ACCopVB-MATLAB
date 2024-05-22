function [U] = sim_posterior_stc_vb(VAMu, VAB, VAD, P, K, R, d, family, LB_dof)
% Sampling one draw of copula data from fitted VB posterior

VAB = reshape(VAB, P, R);
[newParam, ~, ~] = reParam_GaussVB(VAMu, VAB, VAD);
[~, ~, ~, ~, OmegaBar, alpha, delta, nu] = SampPar2STPar(newParam, P, K, d, family, LB_dof);
X = mvstrnd(alpha, OmegaBar, nu, 1);

U = zeros(1,d);
for j = 1:d
    U(:,j) = stcdf(X(:,j),delta(j),nu);
end

end