%%
% logpost = logpost_stc;
recentIndex = floor(iter/10);
f = figure();
plot(logpost_stc(1:recentIndex, 1), logpost_stc(1:recentIndex, 2)) 
title('Plot of skew t copula log-likelihood for every 10 iteration')
xlabel('iteration')
ylabel('log-likelihood')
% exportgraphics(f, 'VB_Likelihood_temp.png')

%%
if iter >= 100
    last_100_VAmu = mean(MatVAMu)';
    [meanG100, ~, ~, ~, meanOmegaBar100, meanAlpha100, meanDelta100, meanNu100] = SampPar2STPar(last_100_VAmu, P, K_G, d, family, LB_dof)
end
clockTime = sum(Time)

%%

addpath("../Empirical_Analysis/")
lambda_hat_Mu = mean(MatVAMu)';
lambda_hat_B = mean(MatVAB)';
lambda_hat_D = mean(MatVAD)';

[OmegaBar_mean_VB, delta_mean_VB, nu_mean_VB, alpha_mean_VB, G_mean_VB, OmegaBar_std_VB, delta_std_VB, nu_std_VB, alpha_std_VB, G_std_VB] = summary_stc_vb(lambda_hat_Mu, lambda_hat_B, lambda_hat_D, P, K_G, K_VAB, d, family, LB_dof);


switch family
    case 'skew-normal'
        nu_mean_VB = 40;
    case 't'
        delta_mean_VB = zeros(size(delta));
    case 'Gaussian'
        nu_mean_VB = 40;
        delta_mean_VB = zeros(size(delta));
end
f = plotquantdepmatrix(delta_mean_VB, reshape(OmegaBar_mean_VB,d,d), nu_mean_VB);
%%
% empirical quantile dependence matrix 
figure()
plot_emp_qtdep(U)

%%
% quantile dependence matrix 
u = 0.05;
quantdepmatrix(u, meanDelta100, meanOmegaBar100, meanNu100)
%%
% Std of VA parameter
% last_100_VAB = reshape(mean(MatVAB(iter-100:iter,:)), P, K_VAB);
% last_100_VAD = diag(mean(MatVAD(iter-100:iter,:)));
% last_100_sqrt = sqrt(diag(last_100_VAB*last_100_VAB' + last_100_VAD.^2));
% last_100_sqrt = [0.3302    0.7255    0.4774    0.2488    0.1414    0.2558    0.9406]';

% [stdG100, ~, ~, ~, stdOmegaBar100, stdAlpha100, stdDelta100, stdNu100] = SampPar2STPar(last_100_sqrt, P, K_G, d)
rng(20220329)
last_100_VAmu = mean(MatVAMu(iter-99:iter,:))';
last_100_VAB = reshape(mean(MatVAB(iter-99:iter,:)), P, K_VAB);
last_100_VAD = mean(MatVAD(iter-99:iter,:))';
sample_T = 5000;

G_sample = zeros(d*K_G,sample_T);
OmegaBar_sample = zeros(d*d,sample_T);
[alpha_sample,delta_sample] = deal(zeros(d, sample_T));
nu_sample = zeros(1, sample_T);
for i = 1:sample_T
    [last_100_Param, ~, ~] = reParam_GaussVB(last_100_VAmu, last_100_VAB, last_100_VAD);
    [G_temp, ~, ~, ~, OmegaBar_temp, alpha_temp, delta_temp, nu_temp] = SampPar2STPar(last_100_Param, P, K_G,  d, family, LB_dof);
    G_sample(:,i) = G_temp(:);
    OmegaBar_sample(:,i) = OmegaBar_temp(:);
    alpha_sample(:,i) = alpha_temp;
    delta_sample(:,i) = delta_temp;
    nu_sample(i) = nu_temp;
end


std_VA = @(mat) sqrt(mean(mat.^2,2) - mean(mat,2).^2);

std_G = reshape(std_VA(G_sample), d, K_G)
std_OmegaBar = reshape(std_VA(OmegaBar_sample), d, d)
% stdG = reshape(sqrt(E_P2(MatG_logdiag) - E2_P(MatG_logdiag)), d, K_G)
std_Alpha = std_VA(alpha_sample)
std_Delta = std_VA(delta_sample)
std_Nu = std_VA(nu_sample)

%%
% windows = 11
meanDelta100 = [   -0.2221   -0.8497   -0.7575]';
meanOmegaBar100 =[
    1.0000   -0.3975   -0.3989
    -0.3975    1.0000    0.8840
    -0.3989    0.8840    1.0000
    ];
meanNu100 = 10.0052;

%%
% windows = 3
meanDelta100 = [   -0.2221   -0.8497   -0.7575]';
meanOmegaBar100 =[
    1.0000   -0.6905   -0.4336
    -0.6905    1.0000    0.3560
    -0.4336    0.3560    1.0000
    ];
meanNu100 = 10.0052;
% plotquantdepmatrix_VIX(meanDelta100, meanOmegaBar100, meanNu100, newX, p);
% plotquantdepmatrix(meanDelta100, meanOmegaBar100, meanNu100)
%%
% windows = 10
meanDelta100 = [   -0.6068   -0.5930   -0.6074]';
meanOmegaBar100 =[
    1.0000   -0.1107   -0.1120
    -0.1107    1.0000    0.8595
    -0.1120    0.8595    1.0000
    ];
meanNu100 = 2;
X_VIX = X(:,1);
step_ahead = 1;

Mdl_VIX = arima(1,0,0);
EstMdl_VIX = estimate(Mdl_VIX,X_VIX,'Display','off');
X_VIX(n+1) = forecast(EstMdl_VIX,step_ahead,X_VIX);

% plotquantdepmatrix_VIX(meanDelta100, meanOmegaBar100, meanNu100, X_VIX, Y(:,1), p, VariableNames);
% plotquantdepmatrix(meanDelta100, meanOmegaBar100, meanNu100, VariableNames)
% [LLMat, URMat, LRMat, ULMat] = quantdepmatrix(0.05, meanDelta100, meanOmegaBar100, meanNu100)
% plot_property('QD_Matrix', 10, 10, 'png')


%%
% Plot Quantile Dependence and print Spearman Correlation
% windows = 10
delta_vb = [
   -0.2953   -0.3752   -0.2625   -0.4101   -0.9405   -0.9488   -0.3101
]';
OmegaBar_vb = [    
    1.0000    0.6437    0.3289    0.4514    0.3646    0.3825    0.4394
    0.6437    1.0000    0.3444    0.4826    0.4480    0.4662    0.4659
    0.3289    0.3444    1.0000    0.4287    0.3798    0.3805    0.3715
    0.4514    0.4826    0.4287    1.0000    0.4946    0.4995    0.4415
    0.3646    0.4480    0.3798    0.4946    1.0000    0.9222    0.4060
    0.3825    0.4662    0.3805    0.4995    0.9222    1.0000    0.4131
    0.4394    0.4659    0.3715    0.4415    0.4060    0.4131    1.0000
];
nu_vb = 3.9518;
plotquantdepmatrix2(delta_vb, OmegaBar_vb, nu_vb, p, VariableNames);
% plot_property('QD_Matrix', 20, 20, 'png')

% [V, Y] = deal(zeros(n,d));
% step_ahead = 1;
% V_hat_1 = zeros(step_ahead,d);
% for j = 1:d
%     Mdl = garch('GARCHLags',1,'ARCHLags',1,'Distribution','t');
%     EstMdl = estimate(Mdl,resid(:,j),'Display','off');
%     [V(:,j), Y(:,j)] = filter(EstMdl,resid(:,j));
%     V_hat_1(:,j) = forecast(EstMdl,step_ahead,'V0',V(:,j));
% %     V1 = infer(EstMdl,resid(:,j));
% end
% V_hat = [V;V_hat_1];
% plotquantdepmatrix3(delta_vb, OmegaBar_vb, nu_vb, p, V_hat, VariableNames);
% [LLMat, URMat, LRMat, ULMat] = quantdepmatrix(0.05, delta_vb, OmegaBar_vb, nu_vb)

% [rhoS] = spearmanrho_AC(delta_vb, OmegaBar_vb, nu_vb)


%%
% Correlation plot
VariableNames = {'.VIX','AAPL.OQ', 'GOOGL.OQ', 'GM.N', 'CAT.N', 'JPM.N', 'BAC.N', 'IBM.N'};

figure
corrplot2(X','Type','Spearman','varNames',VariableNames)
t = title('Skew-t Pseudo Data');
t.FontSize = 15;
figure
corrplot2(U,'Type','Spearman','varNames',VariableNames)
t = title('Copula Data');
t.FontSize = 15;

figure
for i = 1:d
    subplot(d,1,i)
    histogram(X(i,:)','Normalization','probability', 'NumBins', 50)
    tx = xlabel(VariableNames(i));
    ty = ylabel('Density');
    tx.FontSize = 15;
    ty.FontSize = 15;
end
%%
% save('3d_example_04202022_VB', 'MatVAMu', 'MatVAB', 'MatVAD','MatG', 'MatAlpha','MatDelta', 'MatNu', 'MatOmegaBar',...
%     'Time', 'iter')
%%
% For simulation experiment use
% AD: Absolute Difference
% trueAlpha(:)'
% AD_Alpha = abs(trueAlpha(:)' - meanAlpha100)
% % trueDelta(:)'
% AD_Delta = abs(trueDelta(:)' - meanDelta100)
% % trueG(:)'
% AD_G = abs(trueG' - meanG100)
% % trueOmegaBar
% AD_OmegaBar = abs(trueOmegaBar - meanOmegaBar100)
% % trueNu
% AD_Nu = abs(trueNu - meanNu100)