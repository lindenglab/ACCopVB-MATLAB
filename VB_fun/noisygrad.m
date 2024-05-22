function [gradMu, gradB, gradD] = noisygrad(X, parameter, G, V, V1, OmegaBar, alpha, delta, nu, l, w, sigma2, a1, b1, a2, b2, VAMu, VAB, VAD, eps1, eps2)
% Compute noisy gradient of hybrid-ELBO
% INPUTS:
% X: observation data
% G: factor covariance component -- V = G*G' + D
% D: see above
% V: covariance mat
% V1: std component -- OmegaBar = V1*V*V1
% OmegaBar: correlation mat
% alpha: unconstrained skewness parameter vector
% delta: constrained skewness parameter vector
% l: latent variable of skew t distribution -- (0, inf) constrained normal
%   distributed
% w: latent variable of skew t distribution -- Gamma(nu/2, nu/2)
%   distributed
% sigma2: parameters of priors of alpha
% a1, b1: parameters of priors of G
% a2, b2: parameters of priors of nu
% VAMu, VAB, VAD: parameter set lambda of Variational distribution q_lambda
% eps1, eps2: iid r.v from Gaussia Factor VB reparameterization
%
% OUTPUTS:
% gradMu, gradB, gradD: gradients of Gaussian Factor VB parameters lambda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gradLogPost = gradLogPost_TraceGrad01(X, G, V, V1, OmegaBar, alpha, delta, nu, l, w, sigma2, a1, b1, a2, b2);
gradGaussVA = grad_GaussVA(VAMu, VAB, VAD, parameter);

gradMu = gradLogPost - gradGaussVA;
gradB = (gradLogPost - gradGaussVA ) * eps1';
gradB = tril(gradB);
gradD = diag((gradLogPost - gradGaussVA) * eps2');

end

