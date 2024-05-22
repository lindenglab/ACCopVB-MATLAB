function [OmegaBar_mean, delta_mean, nu_mean, alpha_mean, G_mean, OmegaBar_std, delta_std, nu_std, alpha_std, G_std] = summary_stc_vb(VAMu, VAB, VAD, num_var, K, R, d, family, LB_dof, iteration)
% Calculate mean and std of the variational Bayes posterior of the AC
% skew-t model parameters using Monte Carlo sampling
% Note: if family='skew-normal', then nu_mean=100 and nu_std=0
% If family='gaussian', then nu_mean=100, delta=0's, and both std are 0

if nargin < 10
    iteration = 10000;
end
VAB = reshape(VAB, num_var, R);
OmegaBar_mat = zeros(d^2, iteration);
delta_mat = zeros(d, iteration);
nu_mat = zeros(1,iteration);
alpha_mat = zeros(d, iteration);
G_mat = zeros(K*d,iteration);
for j = 1:iteration
    [newParam, ~, ~] = reParam_GaussVB(VAMu, VAB, VAD);
    [G, ~, ~, ~, OmegaBar, alpha, delta, nu] = SampPar2STPar(newParam, num_var, K, d, family, LB_dof);
    OmegaBar_mat(:,j) = OmegaBar(:);
    delta_mat(:,j) = delta;
    nu_mat(j) = nu;
    alpha_mat(:,j) = alpha;
    G_mat(:,j) = G(:);
end
OmegaBar_mean = reshape(mean(OmegaBar_mat,2),d,d);
delta_mean    = mean(delta_mat,2);
nu_mean       = mean(nu_mat);
alpha_mean    = mean(alpha_mat,2);
G_mean        = reshape(mean(G_mat,2),d,K);

fun_std = @(mat) sqrt(mean(mat.^2,2) - mean(mat,2).^2);

OmegaBar_std = reshape(fun_std(OmegaBar_mat),d,d);
delta_std    = fun_std(delta_mat);
nu_std       = fun_std(nu_mat);
alpha_std    = fun_std(alpha_mat);
G_std        = reshape(fun_std(G_mat),d,K);

end