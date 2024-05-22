function [l] = mcmcLatentVarQ(w, x, OmegaBar, delta)
% MCMC Gibbs sampling latent variable q
% 
% INPUT:
% w: latent variable w 
% x: observation data -- d x n matrix
% OmegaBar: correlation matrix of AC's multivariate skew-t distribution
% delta: constrained skewness vector
%
% OUTPUT:
% l: updated l
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps = 1e-6;

Sigma = OmegaBar - delta * delta';
A     =  1 + delta' * Sigma^(-1) * delta;
B     =  w.^(1/2) .* sum(Sigma^(-1)*x.*delta);

l_Mean = B./A;
l_Std  = sqrt(1/A);

l = cnormrnd(l_Mean, l_Std, 0, inf);

counter = 0;
while sum(isinf(l)) > 0
    id = isinf(l);
    l(id) = cnormrnd(l_Mean(id), l_Std, 0, inf);
    counter = counter + 1;
    if counter > 5
        l(id) = eps;
    end
end

end
