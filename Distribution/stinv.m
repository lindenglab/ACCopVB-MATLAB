function [x] = stinv(u, delta, nu, tol)
% Quantile function of AC's skew-t dist
%   u:     nx1, copula data
%   delta: 1,   skewness(constrained)
%   nu:    1,   DoF
%   eps:   1,   precision 

if nargin < 4
    tol = 1e-5;
end

x = zeros(size(u));
% Newton root-finding method
beta = delta2alpha(delta, 1);
for k = 1:300
    Gx = stcdf(x, delta, nu) - u;
    gx = stpdf(x, beta,  nu);
    z = Gx./gx;
    if max(abs(z)) <= tol
        break
    end
%     if max(abs(Gx)) <= eps
%         x = x - z;
%         break
%     end
    x = x - z;
end
%%
% bisection root-finding method
%     fun = @(x) stcdf(x, delta, nu) - u;
%     x = bisection(fun, -10, 10);
end
