function [grad_ThetaTilde_loglike] = grad_st_margin_num(X, ThetaTilde, K_G, family, LB_dof)
% calculalte gradient of cumulative marignal skew-t density (denominator in copula density)

eps = 1e-3;
num_thetatilde = size(ThetaTilde,1);
grad_ThetaTilde_loglike = zeros(num_thetatilde,1);
d = size(X,2);
numG = K_G*d; % # G elements
P = numG + d + 1; % # parameters of skew t dist
for i = num_thetatilde
    ThetaTilde_temp1 = ThetaTilde;
    ThetaTilde_temp2 = ThetaTilde;
    ThetaTilde_temp1(i) = ThetaTilde(i) + eps;
    ThetaTilde_temp2(i) = ThetaTilde(i) - eps;

    [~, ~, ~, ~, ~, ~, delta_temp1, nu_temp1] = SampPar2STPar(ThetaTilde_temp1, P, K_G, d, family, LB_dof);
    [~, ~, ~, ~, ~, ~, delta_temp2, nu_temp2] = SampPar2STPar(ThetaTilde_temp2, P, K_G, d, family, LB_dof);

    loglike_ij_temp1 = loglike_st_margin(X, delta_temp1, nu_temp1);
    loglike_ij_temp2 = loglike_st_margin(X, delta_temp2, nu_temp2);

    loglike_temp1 = sum(loglike_ij_temp1,'all');
    loglike_temp2 = sum(loglike_ij_temp2,'all');
    grad_ThetaTilde_loglike(i) = (loglike_temp1 - loglike_temp2)./ (2*eps);
end

end