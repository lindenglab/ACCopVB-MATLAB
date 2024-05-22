function [x] = stinvspline(u, delta, nu, n)
% if nargin < 4
%     n = 200;
% end
% 
% eps1 = min(u);
% epsN = max(u);
% 
% % q1 = stinv(eps1, delta, nu);
% % qN = stinv(epsN, delta, nu);
% % Q = linspace(q1, qN, n)';
% % P = stcdf(Q, delta, nu);  % skew t cdf based on density function transformation
% % [~, ind] = unique(P); % ind = index of first occurrence of a repeated value 
% % x = interp1(P(ind),Q(ind),u,'spline');
% 
% % % PCHIP method 
% eps2 = eps1+eps1;
% q1 = stinv(eps1, delta, nu);
% qN = stinv(epsN, delta, nu);
% q2 = stinv(eps2, delta, nu);
% Q = linspace(q1, qN, n)';
% % Q = [linspace(q1, q2-eps1, n*0.2), linspace(q2, -q2, n*0.6), linspace(-q2+eps1, qN, n*0.2)]';
% P = stcdf(Q, delta, nu);  % skew t cdf based on density function transformation
% if sum(isinf(Q)) > 0
%     disp("delta: "+delta+", nu: "+nu)
% end
% [~, ind] = unique(P); % ind = index of first occurrence of a repeated value 
% x = interp1(P(ind),Q(ind),u,'pchip');

if nargin < 4
    n = 100;
end
e = 1e-5;

% q1 = stinv(min(u), delta, nu);
% qN = stinv(max(u), delta, nu);
q1 = stinv(e,   delta, nu);
qN = stinv(1-e, delta, nu);

Q = linspace(q1, qN, n)';
% Q = [linspace(q1, 0.199, n*0.4), linspace(0.2, 0.799, n*0.2), linspace(0.8, qN, n*0.4)]';
P = stcdf(Q, delta, nu);  % skew t cdf based on density function transformation

sp = spline(P, Q);
x = ppval(sp, u);

end

