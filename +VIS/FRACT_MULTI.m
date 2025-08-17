function [df, kf, h] = FRACT_MULTI(pp, lbl)
% FRACT_MULTI plots multiple sets of aggregate population data on a...
%   ...fractality domain; i.e., number of primaries vs. normalized...
%   ...gyration diameter.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pp: A cell array of aggregate populations in different conditions,
%       ...each containing primary particle info of individual aggregates.
%   lbl: Labels of different conditions making unique populations of
%       ...aggregates
% ----------------------------------------------------------------------- %
% 
% Outputs:
%   df: fractal dimension (and fitting errors) for each population
%   kf: fractal prefactor ~
%   h: Plot figure handle
% ----------------------------------------------------------------------- %

% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, 500, 750];
set(h, 'color', 'white');

plt = cell(4,1); % initialize plots

ms = [5, 5, 10]; % set marker sizes
mt = {'o', '^', 's'}; % set marker types
% mc = [4,90,141;...
%     116,169,207;
%     208,209,230] / 255; % set of marker colors
% mc = [103,169,207;...
%     239,138,98] / 255;
% mc = [49,163,84;...
%     173,221,142;
%     247,252,185] / 255;
mc = [239,138,98;...
    103,169,207;...
    216,179,101] / 255;

% plot Sorensen's (2011) recommended benchmark
rbar_0 = log10(logspace(0, 1e3, 1e4));
npp_0 = 1.3 * rbar_0.^1.78;
lbl_0 = '$d_\mathrm{f} = 1.78, k_\mathrm{f} = 1.3$';
plt{4} = plot(rbar_0, npp_0, 'Color', [0.4940 0.1840 0.5560], 'LineStyle', '-.',...
    'LineWidth', 4);
hold on

n_set = length(pp); % number of datasets to be plotted

% placeholders for aggregate properties and fit function
npp = cell(n_set,1);
dg = cell(n_set,1);
dpp = cell(n_set,1);
nagg = zeros(n_set,1);
fit = cell(n_set,1);
df = zeros(n_set,3);
kf = zeros(n_set,3);

for i = 1 : n_set
    nagg(i) = length(pp{i}); % number of aggregates in each set    
    
    npp{i} = zeros(nagg(i),1);   
    for j = 1 : nagg(i)
        npp{i}(j) = size(pp{i}{j},1); % number of primaries within each aggregate        
    end
    
    dpp{i} = PAR.GEOMEANPP(pp{i}); % mean primary particle size within aggregates
    dg{i} = PAR.GYRATION(pp{i},npp{i}); % gyration diameter
    
    % plot simulation datapoints
    plt{i} = scatter(dg{i} ./ dpp{i}(:,1), npp{i}, ms(i), mc(i,:), mt{i}, 'LineWidth', 1);
    
    % a fit linear regression weighted with sqrt of number of primaries
    fit{i} = fitlm(table(log(dg{i} ./ dpp{i}(:,1)), log(npp{i})), 'linear',...
        'Weights', sqrt(npp{i})); 
    
    df(i,1) = fit{i}.Coefficients.Estimate(2); % fractal dimension
    kf(i,1) = exp(fit{i}.Coefficients.Estimate(1)); % fractal prefactor
    
    % 95% confidence intervals for fractal properties
    ci = coefCI(fit{i});
    ci_df = ci(2,:);
    ci_kf = exp(ci(1,:));
    
    % 95% ci error bars
    df(i,2) = max(ci_df) - df(1);
    df(i,3) = df(1) - min(ci_df);
    kf(i,2) = max(ci_kf) - kf(1);
    kf(i,3) = kf(1) - min(ci_kf);
        
    % generate the fit data
    npp_fit = kf(i,1) * rbar_0.^df(i,1);
    ci_npp = [ci_kf(1) * rbar_0.^ci_df(1); ci_kf(2) * rbar_0.^ci_df(2)];
    
    % plot the main fit and CI bounds
    plot(rbar_0, npp_fit, 'Color', mc(i,:), 'LineStyle', '--',...
        'LineWidth', 2);
    plot(rbar_0, ci_npp(1,:), rbar_0, ci_npp(2,:), 'Color', mc(i,:),...
        'LineStyle', ':', 'LineWidth', 1);
    
    lbl{i} = strcat(lbl{i}, {' ('}, '$D_\mathrm{f}$ =', {' '}, num2str(df(i,1), '%.2f'),...
        {' '}, '$\pm$', {' '}, num2str(max(df(i,2),df(i,3)), '%.2f'), {','},...
        string(newline), '  $k_\mathrm{f}$ =', {' '}, num2str(kf(i,1), '%.2f'),...
        {' '}, {'$\pm$'}, {' '}, num2str(max(kf(i,2),kf(i,3)), '%.2f'), {')'});
end

box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
legend(cat(2, plt{:})', cat(2, [lbl(:); lbl_0]), 'interpreter', 'latex',...
    'FontSize', 18, 'Location', 'southoutside');
xlabel('$d_\mathrm{g}/d_\mathrm{pp}$ [-]', 'interpreter', 'latex', 'FontSize', 24)
ylabel('$n_\mathrm{pp}$ [-]', 'interpreter', 'latex', 'FontSize', 24)
xlim([2 60])
ylim([4 1200])

% title('Parameteric study of fractality in aggregates', 'FontSize', 22, 'interpreter','latex')
title('(a)', 'FontSize', 24, 'interpreter','latex')

end

