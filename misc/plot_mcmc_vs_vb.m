function [] = plot_mcmc_vs_vb(mcmc_outcome, vb_outcome, titleName, minLimit, maxLimit, x_label, y_label)
% Plot comparsion plot between outcomes from MCMC and VB regarding to mean and standard deviation
mcmcValue = mcmc_outcome(:);
vbValue = vb_outcome(:);
if nargin < 6
    x_label = 'MCMC';
    y_label = 'VB';
end
if nargin < 4
    %     minLimit = min([min(mcmcValue), min(vbValue)]);
    %     maxLimit = max([max(mcmcValue), max(vbValue)]);
%     minLimit = min([mcmcValue,vbValue],[], 'all');
%     maxLimit = max([mcmcValue,vbValue],[], 'all'); 
    minLimit = min([-1,1],[], 'all');
    maxLimit = max([-1,1],[], 'all'); 
end
% figure()

% s = plot(mcmcValue, vbValue, 'diamond', 'MarkerSize',10);
% s.LineWidth = 0.6;
% s.MarkerEdgeColor = 'b';
% s.MarkerFaceColor = [0 0.5 0.5];
s = scatter(mcmcValue, vbValue);
s.Marker = 'o';
s.SizeData = 100;
s.MarkerEdgeColor = [0.4940 0.1840 0.5560]	;
% title_txt = sprintf("%s", titleName);
title_txt = titleName;
t = title(title_txt,'interpreter','latex');
t.FontSize = 20;
xl = xlabel(x_label);
yl = ylabel(y_label);
xl.FontSize = 15;
yl.FontSize = 15;

ylim([minLimit, maxLimit])
xlim([minLimit, maxLimit])

hline = refline(1);
hline.Color = 'k';
axis equal
pbaspect([1 1 1])
end