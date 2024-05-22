function [G, D, V, V1, OmegaBar, alpha, delta, nu] = SampPar2STPar(Param, P, K, d, family, LB_dof)
% Transform sampled paramters (from VB or MCMC) to AC skew-t model parameters
% INPUT:
% Param: input parameter set for re-parameterization
% P: number of parameters
% K: number of copula factors 
% d: number of dimensions
% family: 'skew-t','skew-normal','t','gaussian'
% LB_dof: lower bound of degree of freedom
%
% OUTPUT:
% G, V, V1, D, OmegaBar: factor correlation matrix of AC's skew-t model
%   OmegaBar = V1 * V * V1, V = G*G' + D, V1 = diag(V).^(-1/2)
% alpha: unconstrained skewness vector
% delta: constrained skewness vector
% nu: degrees of freedom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6
    LB_dof = 2;
end
family = lower(family);

diagGIndex = 1:d + 1:d*K; % leading diagonal index of matrix G (required positive values)
LogIndex = [diagGIndex, P]; % logarithmic values index in the parameter vector set
Param(LogIndex) = exp(Param(LogIndex));

G_Num = K * d;
G     = Param(1:G_Num);
alpha = Param(G_Num+1:G_Num+d);
nu    = Param(end) + LB_dof;

G = tril(reshape(G, d, K));
[OmegaBar, D, V, V1] = G2OmegaBar(G);
delta = alpha2delta(alpha, OmegaBar);

switch family
    case 'skew-normal'
        nu = 100;
    case 't'
        [delta, alpha] = deal(zeros(size(delta)));
    case 'gaussian'
        [delta, alpha] = deal(zeros(size(delta)));
        nu = 100;
end

end