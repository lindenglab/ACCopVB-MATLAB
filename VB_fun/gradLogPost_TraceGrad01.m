function [outcome] = gradLogPost_TraceGrad01(X, G, V, V1, OmegaBar, alpha, delta, nu, l, w, sigma2, a1, b1, a2, b2)
% Calculate the gradients of the log-posterior with respect to the parameters listed in Table 5 of the paper.
% The gradients are: grad_logpost_G, grad_logpost_Gtilde, grad_logpost_Alpha, grad_logpost_nu, and grad_logpost_nutilde.
% OUTPUT:
%   gradTG1_5: grad_logpost_G combined with grad_logpost_Gtilde at relative positions
%   gradTalpha1_5: grad_logpost_Alpha
%   gradTnu: grad_logpost_nutilde

n = size(X, 2);
[P_G, K] = size(G);
diagGIndex = 1:P_G + 1:numel(G);

invOmegaBar = inv_OmegaBar(G);
C0 = invOmegaBar + (invOmegaBar*delta*delta'*invOmegaBar)/(1-delta'*invOmegaBar*delta);
C1 = 1 + alpha'*OmegaBar*alpha;
C2 = alpha*alpha'*OmegaBar;
C3 = 1/C1*C2*C0 + 1/C1*C0*C2' - 1/C1^2*C2*C0*C2';
C4 = C0 - C3;
gradTG1 = gradTG(C4, V, V1, G);
gradTalpha1 = gradTAlpha(C0, C1, alpha, OmegaBar);

gradTG2_4 = zeros(size(G));
gradTalpha2_4 = zeros(size(alpha));
for i = 1:n
    C5i = -w(i)*C0*X(:,i)*X(:,i)'*C0;
    C6i = C1^(-1)*C2*C5i + C1^(-1)*C5i*C2' - C1^(-2)*C2*C5i*C2';
    C7i = C5i - C6i;

    C8i = -2*l(i)*w(i)^(1/2)*C1^(-1/2)*C0*X(:,i)*alpha' + l(i)*w(i)^(1/2)*C1^(-3/2)*C2*C0*X(:,i)*alpha';
    C9i = 2*l(i)*w(i)^(1/2)*C0*X(:,i)*delta'*C0;
    C10i = C1^(-1)*C2*C9i + C1^(-1)*C9i*C2' - C1^(-2)*C2*C9i*C2';
    C11i = C8i + C9i - C10i;

    C12i = l(i)^2*C1^(-1/2)*C0*delta*alpha' - 1/2*l(i)^2*C1^(-3/2)*C2*C0*delta*alpha' ...
        + l(i)^2*C1^(-1/2)*alpha*delta'*C0 - 1/2*l(i)^2*C1^(-3/2)*alpha*delta'*C0*C2';
    C13i = - l(i)^2*C0*delta*delta'*C0;
    C14i = C1^(-1)*C2*C13i + C1^(-1)*C13i'*C2' - C1^(-2)*C2*C13i*C2';
    C15i = C12i + C13i - C14i;
    gradTG2_4 = gradTG2_4 + gradTG(C7i+C11i+C15i, V, V1, G);

    C16i = -2*l(i)*w(i)^(1/2)*C1^(-1/2)*OmegaBar*C0*X(:,i) + 2*l(i)*w(i)^(1/2)*C1^(-3/2)*alpha'*OmegaBar*C0*X(:,i)*OmegaBar*alpha;
    C17i = -2*l(i)*w(i)^(1/2)*C0*X(:,i)*delta'*C0;
    C18i = l(i)^2*C1^(-1/2)*OmegaBar*C0*delta - l(i)^2*C1^(-3/2)*OmegaBar*alpha*alpha'*OmegaBar*C0*delta;
    C19i = l(i)^2*C0*delta*delta'*C0;
    gradTalpha2_4 = gradTalpha2_4 + C16i + 2*C18i + gradTAlpha(C5i-C17i-C19i, C1, alpha, OmegaBar);
end

%grad_Gtilde
gradTG5 = -(a1+1) .* sign(G)./(b1+abs(G));
gradTG1_4 = -n/2*gradTG1 - 1/2*gradTG2_4;
gradTG1_5 = gradTG1_4(:) + gradTG5(:);
gradTG1_5(diagGIndex) = gradTG1_5(diagGIndex).*vect(G(diagGIndex)) + 1;

%grad_alpha
gradTalpha5 = -1/sigma2.*alpha;
gradTalpha1_4 = -n/2*gradTalpha1 - 1/2*gradTalpha2_4;
gradTalpha1_5 = gradTalpha1_4(:) + gradTalpha5(:);

%grad_nutilde
gradTnuLogLike = 1/2*sum(log(w)) + n*(1/2*log(nu/2) + 1/2 - 1/2*psi(nu/2))- 1/2*sum(w);
gradTnu = (gradTnuLogLike + (a2-1)/nu - b2)*nu + 1;

outcome = [gradTG1_5; gradTalpha1_5; gradTnu];
end