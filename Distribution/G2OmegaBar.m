function [OmegaBar, D, V, V1] = G2OmegaBar(G)
% MAIN FUNCTION:
% G2OmegaBar: generate multivariate skew t random variables
%
% INPUT:
% mu: mean vector;
% q: q from last iteration
% sigmaW: variance of each w for MH sampling
% x: observation data -- d x n matrix
% OmegaBar: correlation matrix of AC's multivariate skew-t distribution
% delta: constrained skewness vector
% nu: degrees of freedom
%
% OUTPUT:
% w: updated q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 18/08/2021

P_G = size(G,1);
D = speye(P_G);
V = G * G' + D;
V1 = diag(diag(V).^(-1/2));
OmegaBar = V1 * V * V1;

end