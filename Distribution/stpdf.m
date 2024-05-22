function [pdf] = stpdf(x, alpha, nu, islog)
% Compute univariate AC's skew-t density function
% %------------------------------------------------------------------------
% %  Author: Lin Deng, l.deng@mbs.edu
% %========================================================================
% % INPUT:
% % X: 1 x d(dimension) matrix
% % alpha: unconstrained skewness parameter scalar
% % nu: degrees of freedom
% %------------------------------------------------------------------------
% % OUTPUT:
% % pdf: n x 1 density
% %------------------------------------------------------------------------
% % VARIABLE:
% % Q: parallel sandwich matrix multiplication (conditional selection X|Q>0)
% % =======================================================================
% Dec 08, 2023
Q = alpha * x .* sqrt((nu+1)./(nu+x.^2));
if nargin < 4
    pdf = 2 * tpdf(x, nu) .* tcdf(Q, nu+1);
else
    pdf = log(2) + log(tpdf(x, nu)) + log(tcdf(Q, nu+1));
end

end