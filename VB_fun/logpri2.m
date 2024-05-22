function [outcome] = logpri2(G, alpha, nu, sigma2, a1, b1, a2, b2)
    d = size(alpha, 1);
%     nu = nu - 2;
    logG = -(a1+1) * sum(log(1 + abs(G)./b1), 'all');
    logNu = (a2-1)*log(nu) - b2*nu;
    covAlpha = sigma2*eye(d);    
    logAlpha = -1/2.*(alpha'*covAlpha^(-1)*alpha);
    
    outcome = logG + logAlpha + logNu; 
end
