function [w] = mcmcLatentVarW(w, l, sigmaW, X, G, delta, nu)
% RW-MH method approximates latent variable W
% 
% INPUT:
% w: old w from last iteration
% l: latent variable q
% sigmaW: variance of each w for MH sampling
% x: observation data -- d x n matrix
% OmegaBar: correlation matrix of AC's multivariate skew-t distribution
% delta: constrained skewness vector
% nu: degrees of freedom
% 
% OUTPUT:
% w: updated w
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = length(w);
logW     = log(w);
logWProp = normrnd(logW, sqrt(sigmaW));
wProp    = exp(logWProp);

logPostWProp = logpostw(X, G, delta, nu, l, wProp);
logPostW     = logpostw(X, G, delta, nu, l, w);
logMHRatioW  = (logPostWProp - logPostW) + (logWProp - logW);

randList = rand(1, n);
updIndex = randList < min(1, exp(logMHRatioW));

w(updIndex) = wProp(updIndex);
end