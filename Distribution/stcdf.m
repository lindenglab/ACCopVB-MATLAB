function [CDF] = stcdf(x, delta, nu)
%STCDF CDF function of skew-t distribution (A&C)
% INPUT:
% delta: constrained skewness vector
%
% OUTPUT:
% y: CDF
% 24/07/2023 updated

Omega = [1, -delta; -delta, 1];
XU = [zeros(size(x, 1), 1), x];
CDF = 2 * mvtcdf(XU, Omega, nu);
end