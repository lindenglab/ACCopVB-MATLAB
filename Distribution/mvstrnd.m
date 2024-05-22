function [X] = mvstrnd(alpha, OmegaBar, nu, n, option)
% MAIN FUNCTION:
% mvstrnd: generate multivariate skew t random variables
%
% INPUT:
% alpha: unconstrained skewness vectoe
% OmegaBar: correlation matrix
% nu: degrees of freedom
% n: size of desired d dimensional variables (n x d)
%
% OUTPUT:
% X: d dimensional skewT random variable (n x d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 18/08/2021
if nargin < 5
    option = 'B';
else
    option = upper(option);
end

d = size(OmegaBar, 1);
delta = alpha2delta(alpha, OmegaBar);
Omega = cell2mat({OmegaBar, delta;
    delta', 1});

switch option
    case 'A'
        Z_t = mvtrnd(Omega, nu, n);
        L = Z_t(:,end);
        X = Z_t(:,1:d);
        
        idx = L < 0;
        X(idx,:) = - X(idx,:);
    case 'B'
        Z_N = mvnrnd(zeros(1, d+1), Omega, n);
        L = Z_N(:,end);
        X = Z_N(:,1:d);
        
        idx = L < 0;
        X(idx,:) = - X(idx,:);
        
        W = gamrnd(nu/2, 2/nu, [n, 1]);
        X = X.*W.^(-1/2);
end


end