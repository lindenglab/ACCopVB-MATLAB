function [output] = logpostw(x, OmegaBar, delta, nu, q, w)
d = size(x, 1);
Sigma = OmegaBar - delta*delta';
%     output = (d+nu-2)/2.*log(w)...
%               - 1/2.*(w.*sum(Sigma^(-1)*x.*x) - 2*w.^(1/2).*q.*sum(Sigma^(-1)*x.*delta)...
%                + q.^2.*(delta'*Sigma^(-1)*delta) + w.*nu);
output = (d+nu-2)/2.*log(w)...
    - 1/2.*(w.*sum(Sigma^(-1)*x.*x) - 2*w.^(1/2).*q.*sum(Sigma^(-1)*x.*delta)...
    + w.*nu);
end