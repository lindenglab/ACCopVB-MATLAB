function [logpdf, pdf] = stcopula_pdf_implicit(Z, G, alphavec, nu)
% Implicit Copula density of AC's Skew-t
% INPUT:
%   Z:        nxd implicit observations
%   OmegaBar: dxd corr mat
%   deltavec: dx1 skewness(constrained)
%   islog:    0/1 pdf/logpdf
% OUTPUT:
%   PDF:      nx1

[n,d]    = size(Z);
OmegaBar = G2OmegaBar(G);
deltavec = alpha2delta(alphavec, OmegaBar);

P1 = mvstpdf_log(Z, alphavec, OmegaBar, nu, G);
P2 = zeros(n,d);

for j = 1:d
    delta_j = deltavec(j);
    beta_j  = delta2alpha(delta_j, 1);
    P2(:,j) = stpdf_log(Z(:,j), beta_j, nu);
end

logpdf = P1 - sum(P2,2);
pdf    = exp(logpdf);

end