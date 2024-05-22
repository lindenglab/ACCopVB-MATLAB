function [loglike] = loglike_st_margin(X, deltavec, nu)

S        = sqrt((nu+1)./(nu+X.^2));
betavec  = deltavec./sqrt(1-deltavec.^2);
X_margin = betavec'.*S.*X;

P1_margin = gammaln((nu+1)/2) - 1/2*log(nu) - gammaln(nu/2) -(nu+1)/2*log(1+X.^2./nu);
P2_margin = log(tcdf(X_margin, nu+1));

loglike = P1_margin + P2_margin;

end