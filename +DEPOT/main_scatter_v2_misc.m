clc
clear
close all
warning('off')

%% initialize

% parameters for universal correlation
dpp100 = 17.8;
D_TEM = 0.35;

% parameters for Brasil's correlation
alpha_a = 1.08;
k_a = 1.1;

% limits for plotting Brasil's correlation
npp0_lim = [1e0 2e4];
n_npp0 = 1e4;

% limits for plotting universal correlation
da0_lim = [1e0 2e4];
n_da0 = 1e4;

% scattering parameters
mu = 0; % mean
sigma = 0.3; % sd

% Monte Carlo parameters for projected area calculation
n_pnt_mc = 1e2;
n_ang_mc = 5;

%% literature correlations

% initialize scatter data figure
f1 = figure(1);
f1.Position = [50, 50, 1000, 600];
set(f1, 'color', 'white');
tt1 = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile % dpp vs. npp figure

% plot combination of Brasil's and universal correlations
r01 = (npp0_lim(2) / npp0_lim(1)) ^ (1 / (n_npp0 - 1));
npp0 = 1e0 * ones(n_npp0,1) .* r01 .^ (((1 : n_npp0) - 1)');
dpp01 = (((dpp100 ^ (1 / D_TEM)) / 100) * (npp0 / k_a).^(1 / (2 * alpha_a))) .^...
    (D_TEM / (1 - D_TEM));
plt01 = plot(npp0, dpp01, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 2);
hold on

nexttile % dpp vs. da figure

% plot universal correlation
r02 = (da0_lim(2) / da0_lim(1)) ^ (1 / (n_da0 - 1));
da0 = 1e0 * ones(n_da0,1) .* r02 .^ (((1 : n_da0) - 1)');
dpp02 = dpp100 * (da0 / 100) .^ (D_TEM);
plt02 = plot(da0, dpp02, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 2);
hold on

%% non-scattered library data

% load data from first-stage library
load('D:\Hamed\CND\PhD\My Articles\DLCA1\Results\DAT\wsp-1e4.mat', 'pp0')

pars_raw.pp = pp0;
pars_raw.n = zeros(length(pp0), 1);
for i = 1 : length(pp0)
    pars_raw.n(i) = length(pp0{i}); % number of primaries in each aggregate
end
clear pp0

pars_raw = PAR.SIZING(pars_raw); % charcteristic sizes

% projected area diameter
pars_raw.da = 2 * sqrt(PAR.PROJECTION(pars_raw, [], n_pnt_mc, n_ang_mc) / pi);

% plot raw library in dpp vs. npp domain
nexttile(1)
plt1_raw = scatter(pars_raw.n, 1e9 * pars_raw.dpp_g(:,1), 10,...
    hex2rgb('#295F98'), 'o', 'LineWidth', 1);

% plot raw library in dpp vs. da domain
nexttile(2)
plt2_raw = scatter(1e9 * pars_raw.da, 1e9 * pars_raw.dpp_g(:,1), 10,...
    hex2rgb('#295F98'), 'o', 'LineWidth', 1);

%% apply scatter on the perfectly scaled data

% parameters of universal correlation in dpp vs. npp domain, in log-log scale
m = D_TEM / (2 * alpha_a * (1 - D_TEM)); % slope
b = (((dpp100^(1 / D_TEM)) / 100) * ((1 / k_a)^(1 / (2 * alpha_a)))) ^...
    (D_TEM / (1 - D_TEM)); % intercept

% % Compute unit vector perpendicular to the universal dpp vs. npp line
% perp_npp = -1 / sqrt(1 + m^2);  % x-component of the perpendicular direction
% perp_dpp = m / sqrt(1 + m^2);   % y-component of the perpendicular direction

% Generate noise values (dispersion distances) that follow a certain distribution
% dist_scat = normrnd(mu, sigma, size(npp_uc)); % Gaussian
dist_scat = mu + sigma * sqrt(12) * (rand(size(npp_uc)) - 0.5); % uniform
% dist_scat = mu + sigma * randn(size(npp_uc)); % Gaussian

% Scatter the dpp vs. npp points vertically from the universal line in...
% ...dpp vs. npp space
dpp_scat = exp(log(dpp_uc) + dist_scat);
% dpp_scat = exp(log(dpp_uc) + dist_scat / cos(atan(m)));

% rescaling factor for noise around universal correlation
rpp_scat = dpp_scat ./ dpp_uc;

pp_scat = parsdata_sigma{4}(1).pp; % load primary particle data

% introduce noise to the perfectly scales aggregates
for i = 1 : length(pp_scat)
    pp_scat{i}(:,2:5) = repmat(rpp_scat(i), npp_uc(i), 4) .* pp_scat{i}(:,2:5);
end

% make a new aggregate population structure for later calculations
pars_scat.n = npp_uc;
pars_scat.pp = pp_scat;

% recalculate projected area based on the the new scaling
pars_scat.da = 2 * sqrt(PAR.PROJECTION(pars_scat, [], 1e2, 5) / pi);

% calculate GM and GSD of primary particle size
pars_scat.dpp_g = PAR.GEOMEANPP(pars_scat.pp);

%% draw scatterd data

nexttile(1)

plt1_scat = scatter(pars_scat.n, 1e9 * pars_scat.dpp_g(:,1),...
    10, hex2rgb('#EF9C66'), '^', 'LineWidth', 1);

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([1,2000])
ylim([5,50])
xlabel('$n_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)

nexttile(2)

plt2_scat = scatter(1e9 * pars_scat.da, 1e9 * pars_scat.dpp_g(:,1),...
    10, hex2rgb('#EF9C66'), '^', 'LineWidth', 1);

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([10,1500])
ylim([5,50])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)

%% characterize the imposed dispersion in dpp vs. npp domain

fit1 = fitlm(table(log(pars_scat.n), log(1e9 * pars_scat.dpp_g(:,1))),...
    'linear', 'Weights', sqrt(pars_scat.n)); % fit a linear regression...
% ...weighted by sqrt of number of primaries

D_beta_scat = fit1.Coefficients.Estimate(2); % ensemble scaling exponent
k_beta_scat = exp(fit1.Coefficients.Estimate(1)); % ensemble scaling prefactor

% 95% confidence intervals for fractal properties
ci1 = coefCI(fit1);
ci_D_beta_scat = ci1(2,:);
ci_k_beta_scat = exp(ci1(1,:));

% 95% ci error bars
dci_D_beta_scat = max(ci_D_beta_scat) - D_beta_scat;
dcip_k_beta_scat = max(ci_k_beta_scat) - k_beta_scat;
dcin_k_beta_scat = k_beta_scat - min(ci_k_beta_scat);

% generate the fit data
dpp_scat_fit1 = k_beta_scat * (npp0.^D_beta_scat);
ci_dpp_scat_fit1 = [ci_k_beta_scat(1) * (npp0.^ci_D_beta_scat(1)),...
    ci_k_beta_scat(2) * (npp0.^ci_D_beta_scat(2))];

nexttile(1)
% plot the main fit and CI bounds
plt1_fit = plot(npp0, dpp_scat_fit1, 'Color', hex2rgb('#C96868'),...
    'LineStyle', '--', 'LineWidth', 1.5);
plt1_err = plot(npp0, ci_dpp_scat_fit1(:,1), npp0, ci_dpp_scat_fit1(:,2),...
    'Color', hex2rgb('#C96868'), 'LineStyle', ':', 'LineWidth', 1);

lgd_uc1 = strcat(string(newline), 'Olfert $\&$ Rogak (2019)', string(newline),...
    'and Brasil et al. (1991)', string(newline), '$D_\mathrm{\beta_{uc}}$ =',...
    {' '}, num2str(m, '%.2f'), ', $k_\mathrm{\beta_{uc}}$ =', {' '}, num2str(b, '%.2f'));

lgd_fit1 = strcat(string(newline), string(newline),...
    '$D_\mathrm{\beta_{scat}}$ =', {' '}, num2str(D_beta_scat, '%.2f'),...
    {' '}, '$\pm$', {' '}, num2str(dci_D_beta_scat, '%.2f'), {','},...
    string(newline), '$k_\mathrm{\beta_{scat}}$ =', {' '}, num2str(k_beta_scat, '%.2f'),...
    {' '}, {'$\pm$'}, {' '}, num2str(max(dcip_k_beta_scat, dcin_k_beta_scat), '%.2f'));

legend(cat(2, plt1_raw, plt1_scat, plt01, plt1_fit),...
    cat(2, {'Original scaling'}, {'Secondary dispersion'},...
    lgd_uc1, lgd_fit1), 'interpreter', 'latex', 'FontSize', 11,...
    'location', 'southoutside', 'orientation', 'horizontal', 'NumColumns', 2)


%% check scaling properties for dpp vs. da

fit2 = fitlm(table(log(1e7 * pars_scat.da), log(1e9 * pars_scat.dpp_g(:,1))),...
    'linear', 'Weights', sqrt(pars_scat.n)); % fit a linear regression...
% ...weighted by sqrt of number of primaries

D_TEM_scat = fit2.Coefficients.Estimate(2); % ensemble scaling exponent
dpp100_scat = exp(fit2.Coefficients.Estimate(1)); % ensemble scaling prefactor

% 95% confidence intervals for fractal properties
ci2 = coefCI(fit2);
ci_D_TEM_scat = ci2(2,:);
ci_dpp100_scat = exp(ci2(1,:));

% 95% ci error bars
dci_D_TEM_scat = max(ci_D_TEM_scat) - D_TEM_scat;
dcip_dpp100_scat = max(ci_dpp100_scat) - dpp100_scat;
dcin_dpp100_scat = dpp100_scat - min(ci_dpp100_scat);

% generate the fit data
dpp_scat_fit2 = dpp100_scat * ((da0/100).^D_TEM_scat);
ci_dpp_scat_fit2 = [ci_dpp100_scat(1) * (da0/100).^ci_D_TEM_scat(1),...
    ci_dpp100_scat(2) * (da0/100).^ci_D_TEM_scat(2)];

nexttile(2)
% plot the main fit and CI bounds
plt2_fit = plot(da0, dpp_scat_fit2, 'Color', hex2rgb('#C96868'),...
    'LineStyle', '--', 'LineWidth', 1.5);
plt2_err = plot(da0, ci_dpp_scat_fit2(:,1), da0, ci_dpp_scat_fit2(:,2),...
    'Color', hex2rgb('#C96868'), 'LineStyle', ':', 'LineWidth', 1);

lgd_uc2 = strcat(string(newline), string(newline), 'Olfert $\&$ Rogak (2019)',...
    string(newline), '$D_\mathrm{TEM}$ = 0.35, $d_\mathrm{pp,100}$ = 17.8');

lgd_fit2 = strcat(string(newline), string(newline), '$D_\mathrm{TEM_{scat}}$ =',...
    {' '}, num2str(D_TEM_scat, '%.2f'), {' '}, '$\pm$', {' '},...
    num2str(dci_D_TEM_scat, '%.2f'), {','}, string(newline),...
    '$d_\mathrm{pp,100_{scat}}$ =', {' '}, num2str(dpp100_scat, '%.2f'),...
    {' '}, {'$\pm$'}, {' '}, num2str(max(dcip_dpp100_scat, dcin_dpp100_scat), '%.2f'));

legend(cat(2, plt2_raw, plt2_scat, plt02, plt2_fit),...
    cat(2, {'Original scaling'}, {'Secondary dispersion'},...
    lgd_uc2, lgd_fit2), 'interpreter', 'latex', 'FontSize', 11,...
    'location', 'southoutside', 'orientation', 'horizontal', 'NumColumns', 2)

%% Draw resulting aggregates

f2 = figure(2);
f2.Position = [100, 100, 1200, 800];
set(f2, 'color', 'white');

tt2 = tiledlayout(2, 3, 'Padding', 'none', 'TileSpacing', 'none');

jj_render = sort(randperm(length(pp_scat), 6));

for j = 1 : length(jj_render)

    % figure(j+2)
    % set(gcf, 'Position', [(2*j+3)*50, (2*j-1)*50, 500, 500], 'color', 'white');

    nexttile
    UTILS.PLOTPP(pp_scat{j}(:,3), pp_scat{j}(:,4), pp_scat{j}(:,5),...
        pp_scat{j}(:,2))

end

%% Export plots

exportgraphics(f1, 'outputs\dispersion.png',...
    'BackgroundColor','none', 'ContentType','vector', 'Resolution', 150)
exportgraphics(f2, 'outputs\render.png',...
    'BackgroundColor','none', 'ContentType','vector', 'Resolution', 150)
