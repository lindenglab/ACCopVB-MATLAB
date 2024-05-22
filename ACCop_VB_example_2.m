%%
clc;clear;

% Generate random number and add the relevant directory paths
rng(20240101)
addpath('VB_fun/')
addpath('misc/')
addpath('Data/')
addpath('Distribution/')

% Loading AUD Exchange Rate data - this dataset is publically available and is downloaded from, 
% https://www.rba.gov.au/statistics/historical-data.html#exchange-rates
X = readtable("Data/2018-2022-num.xls");
X = X{:,:};
X_return = diff(log(X));
[n,d] = size(X_return);

% Converting copula data by using kernel cdf
% The distribution object is saved for construction of predictive
% distributions.
U = zeros(n,d);
pd_cell = {d};
for j = 1:d
    pd = fitdist(X_return(:,j),'Kernel');
    U(:,j) = cdf(pd, X_return(:,j));
    pd_cell{j} = pd;
end

% parsave('Data/example_2_kernel',pd_cell,U);
%% VB Estimation
% Settings for the VB algorithm here
% K=number of Gaussian factors in the covariance matrix of the AC skew-t copula model (see equation xx)
% R=number of factors in the covariance matrix of the Gaussian variational approximation
K = 2;
R = 3;

% T=number of optimization steps for the VB approximation
% M=number of MCMC iteration to sample the latent variables within each VB optimization step
T = 5000;
M = 20;

% family is the type of copula to estimate: 'skew-t' (default),
% 'skew-normal' or 'gaussian'. The 'skew-normal' is the implicit copula
% constructed from AC skew-normal distribution
% LB_dof is the lower bound for the degree of freedom (nu) for the copula
family = 'skew-t';
LB_dof = 2;

% VB estimation of the AC skew-t copula
[Lambda,logpost_stc,Time] = vb_st_copula_opt_b(U, K, R, T, M, 1, family, LB_dof);

% Summarizing the VB fit 
VAMu = Lambda.VAMu;
VAB  = Lambda.VAB;
VAD  = Lambda.VAD;

d = size(U,2);
P = d+d*K+1;

% Calculating the mean and standard deviation of the model parameters:
% Omega, delta, nu, alpha and G (notation as per the AC skew-t defined in
% the paper, Section 2)
[OmegaBar_mean, delta_mean, nu_mean, alpha_mean, G_mean, OmegaBar_std, delta_std, nu_std, alpha_std, G_std] = summary_stc_vb(VAMu, VAB, VAD, P, K, R, d, family, LB_dof);

%% Saving the results using the current date
datetxt = datetime("now",'Format','ddMMMyyyy_HHmmss');
destinationFolder = sprintf('Results_%s', datetxt)
if ~exist(destinationFolder, 'dir')
    mkdir(destinationFolder);
end

filename = sprintf('%s/d%d_K%d_R%d.mat', destinationFolder, d, K, R);
parsave(filename,family,LB_dof,d,P,K,R,Lambda,logpost_stc,Time,...
    OmegaBar_mean, delta_mean, nu_mean, alpha_mean, G_mean, OmegaBar_std, delta_std, nu_std, alpha_std, G_std);

%% Plotting summary of quantile dependence metrics from saved results
% Note that the data and VB file can be changed
addpath('Results/'); %directory containing results from paper

% Loading the VB fit
Fit_File = load("d17_K2_R3.mat");

% Calculating the empirical and analytical skew-t quantile dependence metrics
% with quantile level at u 
u = 0.05;
[LL_Ori, UR_Ori, LR_Ori, UL_Ori] = empQuantDepMat(U, u);
[LL_Fit, UR_Fit, LR_Fit, UL_Fit] = stQuantDepMat(u,Fit_File.delta_mean, Fit_File.OmegaBar_mean, Fit_File.nu_mean);

% Generating figure for comparing quantile dependence
figure()
subplot(2,2,1)
plot_mcmc_vs_vb(LL_Ori,LL_Fit,'$\lambda_{LL}$',0,1,"True","VB")
subplot(2,2,2)
plot_mcmc_vs_vb(UR_Ori,UR_Fit,'$\lambda_{UR}$',0,1,"True","VB")
subplot(2,2,3)
plot_mcmc_vs_vb(LR_Ori,LR_Ori,'$\lambda_{LR}$',0,1,"True","VB")
subplot(2,2,4)
plot_mcmc_vs_vb(UL_Ori,UL_Ori,'$\lambda_{UL}$',0,1,"True","VB")

%% Sample from the posterior predictive distribution of the copula, with uncertainty in model unknowns integrated out
%This can be used to obtain predictive distributions of the base variable
addpath('Results/')
Fit_File = load("d17_K2_R3.mat");

P = Fit_File.P;
K = Fit_File.K;
R = Fit_File.R;
d = Fit_File.d;

Lambda = Fit_File.Lambda;
family = Fit_File.family;
LB_dof = Fit_File.LB_dof;

%Simulating 10k draws from the posterior predictive
Iteration = 1e4;
U_mat = zeros(Iteration,d);
X_mat = zeros(Iteration,d);
pd_cell = load('Data/example_2_kernel.mat').pd_cell;

for i = 1:Iteration
    U_sample = sim_posterior_stc_vb(Lambda.VAMu, Lambda.VAB, Lambda.VAD, P, K, R, d, family, LB_dof);
    U_mat(i,:) = U_sample;

    %% Transform the simulated copula data to the raw data scale by inverting the
    % % specified marginal distribution. eg: kernel density
    % % Marginal models can also be a dynamic univariate model
    X_sample = zeros(size(U_sample));
    for j = 1:d
        X_sample(:,j) = icdf(pd_cell{j},U_sample(:,j));
    end
    X_mat(i,:) = X_sample;
end

