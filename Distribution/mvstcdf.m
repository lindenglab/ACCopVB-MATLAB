function [CDF] = mvstcdf(X, delta, OmegaBar, nu)
%     Omega = [1, rho; rho, 1]; % because of the order of Omega, we still
%     follow the order X1 = [0, X11, X12]
% Omega = {1, -delta';
%     -delta, OmegaBar};
% Omega = cell2mat(Omega);
% XU = [zeros(size(X, 1), 1), X];
% Y = 2 * mvtcdf(XU, Omega, nu);

n = size(X,1);
Omega = {OmegaBar, -delta;
        -delta',   1};
Omega = cell2mat(Omega);
X2 = [ X, zeros(n,1) ];
CDF = 2 * mvtcdf(X2, Omega, nu);

end