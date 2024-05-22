function [log_PDF] = mvstpdf_log(X, alpha, OmegaBar, nu, G)
% Compute log density of multivariate AC's skew-t distribution
% %------------------------------------------------------------------------
% %  Author: Lin Deng, l.deng@mbs.edu
% %========================================================================
% % INPUT:
% % X: n x d(dimension) matrix
% % OmegaBar: correlation matrix
% % alpha: unconstrained skewness parameter vector
% % nu: degrees of freedom
% %------------------------------------------------------------------------
% % OUTPUT:
% % PDF: n x 1 density
% %------------------------------------------------------------------------
% % VARIABLE:
% % Q: parallel sandwich matrix multiplication (conditional selection X|Q>0)
% % P: input n x 1 matrix for univariate t CDF
% % =======================================================================
% Jan 16, 2023

if nargin < 5
    invOmegaBar = OmegaBar^(-1);
else
    invOmegaBar = inv_OmegaBar(G);
end

d = size(X, 2);
Q = sum(invOmegaBar * X' .* X', 1);
P = alpha' * X' .* sqrt((nu + d) ./ (nu + Q));


P1 = log(gamma((nu+d)/2)) - log(gamma(nu/2)) - d/2*log(nu) - 1/2*log(det(OmegaBar)) - d/2*log(pi);
P2 = -(nu+d)/2 * ( 1 + 1/nu * sum( invOmegaBar * X' .* X') );
P3 = log( tcdf(P', nu + d) );
log_PDF = P1 + P2' + P3;
% PDF = 2 * mvtpdf(X, OmegaBar, nu) .* tcdf(P', nu + d);

end