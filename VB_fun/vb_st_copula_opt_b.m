function [Lambda, logpost_stc, Time] = vb_st_copula_opt_b(U, K, R, T, M, Verbose, family, LB_dof, init_Lambda)
% Variational Bayes Method for Skew T copula Model
% INPUT
% U: copula data
% K: Gaussian factor AC-skew t covariance is d-by-K
% R: Gaussian factor variational covariance is P-by-R
% T: num of iterations for VB process
% M: num of iterations for latent variables' MCMC processes
% Verbose: 1 indicates return converging time for every 100 iteration
% family: family selection includes "skew-t", "skew-normal", "t", "gaussian"
% LB_dof: lower bound of degree of freedom
% init_Lambda: initial value of VB parameters
%
% OUTPUT
% Lambda: fitted VB parameters
% logliks_stcopula: log-likelihood of selected copula model
% Time: running time for VB processes

[n, d] = size(U);
X      = zeros(d, n);
Time   = 0; 

numG = K*d;        % # G elements
P    = numG + d + 1; % # parameters of skew t dist

if nargin < 7
    family = 'skew-t';
end
family = lower(family);

if nargin < 8
    LB_dof = 2;
end

if nargin < 9
    nu = 25;
    VAMu      = zeros(P,1);
    VAMu(end) = log(nu);
    VAB       = 1e-3*tril(ones(P,R));
    VAD       = 1e-3*ones(P,1);
else
    VAMu = init_Lambda.VAMu;
    VAB  = reshape(init_Lambda.VAB, P, R);
    VAD  = init_Lambda.VAD;
    [newParam, ~, ~] = reParam_GaussVB(VAMu, VAB, VAD);
    [~, ~, ~, ~, ~, ~, ~, nu] = SampPar2STPar(newParam, P, K, d, family, LB_dof);
end

sigmaW = 1*ones(1, n); 
w_vec  = gamrnd(nu/2, 2/nu, 1, n);
l_vec  = cnormrnd(zeros(1,n), ones(1,n), 0, inf);

% Prior dist parameter
sigma2 = 25; % alpha ~ N(0, sigma2)
a1 = 3; % G ~ GDP(3, 1)
b1 = 1;
a2 = 3; % Nu ~ f_gamma(3, rate = 0.2)
b2 = 0.2;

% StoreRoom
[stateMu, stateB, stateD] = deal(struct);
[MatVAD, MatVAMu] = deal(zeros(500, P));
MatVAB = zeros(500, R*P);
logpost_stc = zeros(round(T/10),2);

%%
for iter = 1:T
    tic

    [newParam, eps1, eps2] = reParam_GaussVB(VAMu, VAB, VAD);
    [G, ~, V, V1, OmegaBar, alpha, delta, nu] = SampPar2STPar(newParam, P, K, d, family, LB_dof);

    for j = 1:d
        X(j, :) = stinvspline(U(:, j), delta(j), nu);
    end


    for i = 1:M
        w_vec = mcmcLatentVarW(w_vec, l_vec, sigmaW, X, OmegaBar, delta, nu);
        l_vec = mcmcLatentVarQ(w_vec, X, OmegaBar, delta);
    end

    [gradMu, gradB, gradD] = noisygrad(X, newParam, G, V, V1, OmegaBar, alpha, delta, nu, l_vec, w_vec, sigma2, a1, b1, a2, b2, VAMu, VAB, VAD, eps1, eps2);
    
    [updMu, stateMu] = Adadelta(gradMu, stateMu);
    [updB,  stateB]  = Adadelta(gradB, stateB);
    [updD,  stateD]  = Adadelta(gradD, stateD);
    
    VAMu = VAMu + updMu;
    VAB  = tril(VAB + updB);
    VAD  = VAD + updD;

    if ~mod(iter, 10)
        logLike = sum(stcopula_pdf_implicit(X', G, alpha, nu));
        logPri  = logpri2(G, alpha, nu, sigma2, a1, b1, a2, b2);

        logpost_stc(iter/10,1) = iter;
        logpost_stc(iter/10,2) = logLike + logPri;
    end

    Time = Time + toc;
    if ~mod(iter,100) && Verbose
        avgTime = Time/iter;
        disptxt = sprintf("VB Iteration: %d/%d || AvgTime: %.4fs || TotalTime: %.4fmin \n", iter, T, avgTime, Time/60);
        fprintf(disptxt)
    end

    if iter <= 500
        MatVAMu(iter, :) = VAMu(:)';
        MatVAB(iter, :)  = VAB(:)';
        MatVAD(iter, :)  = VAD(:)';

    else
        MatVAMu(1:end-1,:) = MatVAMu(2:end,:);
        MatVAB(1:end-1,:)  = MatVAB(2:end,:);
        MatVAD(1:end-1,:)  = MatVAD(2:end,:);

        MatVAMu(end, :) = VAMu(:)';
        MatVAB(end, :)  = VAB(:)';
        MatVAD(end, :)  = VAD(:)';
    end

end

Lambda = struct;
if T <= 500
    Lambda.VAMu = VAMu(:);
    Lambda.VAB  = VAB(:);
    Lambda.VAD  = VAD(:);
else
    Lambda.VAMu = mean(MatVAMu)';
    Lambda.VAB  = mean(MatVAB)';
    Lambda.VAD  = mean(MatVAD)';
end

end