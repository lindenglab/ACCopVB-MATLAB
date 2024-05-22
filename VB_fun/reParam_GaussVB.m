function [parameter, eps1, eps2] = reParam_GaussVB(VAMu, VAB, VAD)
% Reparameterize Gaussian Factor Variational Bayes parameters
%
% INPUT:
% VAMu: Gaussian VA mu
% VAB: Gaussian VA B
% VAD: Gaussian VA D (diagonal value vector)
% 
% OUTPUT:
% parameter: derived from re-parameterization
% eps1: VAB_K x 1 vector sampled from N(0, 1)
% eps2: VAB_P x 1 vector sampled from N(0, 1), VAB_P = P: #
%       distribution parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[P, VAB_K] = size(VAB);
eps1 = normrnd(0, 1, [VAB_K, 1]);
eps2 = normrnd(0, 1, [P, 1]);

parameter = VAMu + VAB*eps1 + VAD.*eps2;
end