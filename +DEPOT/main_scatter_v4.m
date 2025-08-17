clc
clear
close all
warning('off')

%% Initialize and define key parameters %%

% Olfert and Rogak (2019)'s universal correlation
D_TEM = 0.35; % exponent
dpp100 = 17.8; % pefactor
da_lim_uc = [1e0 2e4];  % limits on the projected area diameter (for plotting)
n_da_uc = 1e4; % plot data counts
uc = @(y) dpp100 * (y / 100) .^ D_TEM; % on-demand function for the...
    % ...forward correlation (da and dpp are in [nm])
uc_inv = @(x) 100 * (x / dpp100) .^ (1 / D_TEM); % inverse function to get...
    % ...da from dpp
b0_uc = dpp100 * 100^(-D_TEM);

% Brasil et al. (1999)'s correlation
alpha_a = 1.08; % exponent
k_a = 1.1; % pefactor
npp_lim_bc = [1e-1 2e4]; % limits on the number of primaries (for plotting)
n_npp_uc = 1e4; % plot data counts
bc = @(x,y) k_a * (x./y).^(2 * alpha_a); % on-demand function for the...
    % ...forward correlation
bc_inv = @(z) (z / k_a).^(1 / (2 * alpha_a)); % inverse function to get...
    % ...dpp/da from npp

% combination of universal and Brasil correlations
m_buc = D_TEM / (2 * alpha_a * (1 - D_TEM)); % exponent
b0_buc = (((dpp100^(1 / D_TEM)) / 100) * ((1 / k_a)^(1 / (2 * alpha_a)))) ^...
    (D_TEM / (1 - D_TEM)); % prefactor
buc = @(z) b0_buc * (z .^ m_buc); % convert, on-demand, number of...
    % ...primaries to primary particle size based on the correlations 
buc_inv = @(x) (b0_buc * x) .^ (1 / m_buc); % inverse combined...
    % ...function to get npp from dpp

% spread variables for the baseline aggregate size distribution
gm_da = 100; % geometric mean of projected area diameter distribution
gsd_da = 1.6; % geometric standard deviation of ~
% mu_da = log(gm_da^2 / sqrt(gsd_da + gm_da^2)); % convert to mean in log-log space
% sigma_da = sqrt(log(1 + (gsd_da / gm_da^2))); % sd in log-log space

% spread variables for the noise
mu_scat = 0; % Gaussiam mean of noise distribution around the universal...
    % ...correlation in log-log space
sigma_scat = 0.2; % Gaussian standard deviation
da_lim_noise = [2e1 2e3]; % extents of noise generation 
cn_scat = 0.6;

% address of aggregate library to be imported for scaling and dispersion
fdir = 'D:\Hamed\CND\PhD\My Articles\DLCA1\Results\DAT\Desktop-simulations\AUG-02-22\sigmapp1\3';
fname = 'wsp_sigmapp_1.3_Aug9';
varname = 'pp0';
vardir = '';

% resolution for projected area calculation
n_mc = 1e4;
n_ang = 20;

%% Raw data against the Brasil's and universal correlations %%

% load non-scaled first-stage library data
load(strcat(fdir, '\', fname, '.mat'), varname)

% create particle structure  
pars_raw.pp = eval(strcat(varname, vardir)); % store primary particle info
n_agg_raw = length(pars_raw.pp); % number of aggregates
pars_raw.n = zeros(n_agg_raw,1);
for i = 1 : n_agg_raw
    pars_raw.n(i) = size(pars_raw.pp{i}, 1); % count the number of primaries
end
pars_raw = PAR.SIZING(pars_raw); % get characteristic sizes

eval(['clear ', varname])

% initialize figure for raw data vs. correlations
f1 = figure(1);
f1.Position = [50, 50, 800, 450];
set(f1, 'color', 'white');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact')

nexttile(1) % dpp vs. npp figure

% plot combination of Brasil's and universal correlations
r_bc = (npp_lim_bc(2) / npp_lim_bc(1)) ^ (1 / (n_npp_uc - 1));
npp_bc = npp_lim_bc(1) * ones(n_npp_uc,1) .* r_bc .^ (((1 : n_npp_uc) - 1)');
dpp_bc = buc(npp_bc);
% dpp_bc = (((dpp100 ^ (1 / D_TEM)) / 100) * (npp_bc / k_a).^(1 / (2 * alpha_a))) .^...
%     (D_TEM / (1 - D_TEM));
plt1a_bc = plot(npp_bc, dpp_bc, 'Color', hex2rgb('#597445'),...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

% plot stage-1 library of aggregates in dpp vs npp domain
plt1a_raw = scatter(pars_raw.n, 1e9 * pars_raw.dpp_g(:,1), 8,...
    hex2rgb('#CDC2A5'), 'o', 'LineWidth', 1);

% appearance configs for dpp vs npp subplot
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([0.8 * min(pars_raw.n), 1.2 * max(pars_raw.n)])
ylim([1e9 * 0.9 * min(pars_raw.dpp_g(:,1)), 1e9 * 1.1 * max(pars_raw.dpp_g(:,1))])
xlabel('$n_\mathrm{pp}$ [-]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
box on

nexttile(2) % dpp vs. da figure

% plot universal correlation
r_uc = (da_lim_uc(2) / da_lim_uc(1)) ^ (1 / (n_da_uc - 1));
da_uc = da_lim_uc(1) * ones(n_da_uc,1) .* r_uc .^ (((1 : n_da_uc) - 1)');
dpp_uc = uc(da_uc);
plt1b_uc = plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

% plot library in dpp vs. da domain
pars_raw.da = 2 * sqrt(PAR.PROJECTION(pars_raw, [], n_mc, n_ang) / pi);
plt1b_raw = scatter(1e9 * pars_raw.da, 1e9 * pars_raw.dpp_g(:,1), 8,...
    hex2rgb('#CDC2A5'), 'o', 'LineWidth', 1);

% appearance configs for dpp vs da subplot
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([1e9 * 0.9 * min(pars_raw.da(:,1)), 1e9 * 1.1 * max(pars_raw.da(:,1))])
ylim([1e9 * 0.9 * min(pars_raw.dpp_g(:,1)), 1e9 * 1.1 * max(pars_raw.dpp_g(:,1))])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
box on

% generate legend for the entire tile
lgd1 = legend(cat(2, plt1a_bc, plt1b_uc, plt1a_raw),...
    cat(2, {strcat('Brasil et al. (1999) +', string(newline),...
    'Olfert $\&$ Rogak (2019)')}, {'Olfert $\&$ Rogak (2019)'},...
    {strcat('First-stage simulated aggregates', string(newline),...
    '(non-scaled, non-filtered)')}), 'interpreter', 'latex', 'FontSize', 11,...
    'orientation', 'horizontal', 'NumColumns', 3);
lgd1.Layout.Tile = 'south';

%% Generate Gaussian random noise around the universal correlation %%

n_scat = round(cn_scat * length(pars_raw.n)); % number of noise points for random selection

% % uniform baseline distribution of projected area size for scattering to be...
%     % ...implemented
% r_scat = (da_lim_noise(2) / da_lim_noise(1)) ^ (1 / (n0_scat - 1));
% da0_scat = da_lim_noise(1) * ones(n0_scat,1) .* r_scat .^ (((1 : n0_scat) - 1)');

% make a log-normal distribution of projected area diameter
% da0_scat = lognrnd(mu_da, sigma_da, [n0_scat, 1]);
da_scat = exp(normrnd(log(gm_da), log(gsd_da), [n_scat, 1]));

dpp00_scat = uc(da_scat); % convert to baseline primary particle size

% Compute unit vector perpendicular to the universal dpp vs. da line
% h2_perp_da = -1 / sqrt(1 + D_TEM^2);  % x-component
% h2_perp_dpp = D_TEM / sqrt(1 + D_TEM^2);  % y-component

% Generate noise values (dispersion distances) that follow a certain distribution
% dist_scat = zeros(n0_scat, 1); % test out with zero noise
% dist_scat = mu + sigma * sqrt(12) * (rand(size(n_noise)) - 0.5); % uniform
dist_scat = mu_scat + sigma_scat * randn(n_scat, 1); % Gaussian
% dist_scat = normrnd(mu_scat, sigma_scat, [n0_scat, 1]); % Gaussian

% Scatter the points across the universal dpp vs. da correlation
% da_scat = exp(log(da0_scat) + dist_scat * h2_perp_da);
% dpp0_scat = exp(log(dpp00_scat) + dist_scat * h2_perp_dpp);
dpp0_scat = exp(log(dpp00_scat) + dist_scat);

npp0_scat = bc(da_scat, dpp0_scat); % raw converted number of primaries

npp_scat = round(npp0_scat); % corrected number of primaries (has to be integer)

dpp_scat = da_scat ./ bc_inv(npp_scat); % get corrected primary particle size

% initialize figure for baseline monte carlo points
f2 = figure(2);
f2.Position = [100, 100, 900, 600];
set(f2, 'color', 'white');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact')

nexttile(1) % dpp vs. npp figure

% plot combination of Brasil's and universal correlations
plt2a_bc = plot(npp_bc, dpp_bc, 'Color', hex2rgb('#597445'),...
    'LineStyle', '-.', 'LineWidth', 3);
lgd2a_bc = strcat('Brasil et al. (1999) + Olfert $\&$ Rogak (2019)',...
    string(newline), '$D_\mathrm{BUC}$ =', {' '}, num2str(m_buc, '%.2f'),...
    ', $k_\mathrm{BUC}$ =', {' '}, num2str(b0_buc, '%.2f'));
hold on

% plot random gussian points generated above in dpp vs npp domain
plt2a_scat = scatter(npp0_scat, dpp0_scat, 12, hex2rgb('#789DBC'), '.',...
    'LineWidth', 1);
lgd2a_scat = strcat(string(newline), 'Transformation of Gaussian noise seeds', string(newline),...
    'from $d_\mathrm{pp}$-$d_\mathrm{a}$ space to $d_\mathrm{pp}$-$n_\mathrm{pp}$ space',...
    string(newline), 'using correlation of Brasil et al. (1999)');

% appearance configs for dpp vs npp subplot
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([0.8 * min(npp0_scat), 1.2 * max(npp0_scat)])
ylim([0.9 * min(dpp0_scat), 1.1 * max(dpp0_scat)])
xlabel('$n_\mathrm{pp}$ [-]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
box on

nexttile(2) % dpp vs. da figure

% plot universal correlation
plt2b_uc = plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
lgd2b_uc = strcat('Olfert $\&$ Rogak (2019)',...
    string(newline), '$D_\mathrm{TEM}$ =', {' '}, num2str(D_TEM, '%.2f'),...
    ', $k_\mathrm{TEM}$ =', {' '}, num2str(b0_uc, '%.2f'));
hold on

% plot library in dpp vs. da domain
plt2b_scat = scatter(da_scat, dpp0_scat, 12, hex2rgb('#B06161'), '.',...
    'LineWidth', 1);
lgd2b_scat = strcat(string(newline), 'Gaussian noise seeds (in log-log space)',...
    string(newline), '$\mu$ =', {' '}, num2str(mu_scat, '%.0f'), ',', {' '},...
    '$\sigma$ =', {' '}, num2str(sigma_scat, '%.1f'));

% appearance configs for dpp vs da subplot
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([0.9 * min(da_scat), 1.1 * max(da_scat)])
ylim([0.9 * min(dpp0_scat), 1.1 * max(dpp0_scat)])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
box on

%% evaluate trends of "non-rounded" Monte Carlo points  %%

figure(f2);

%%% find the best curve fit to the Monte Carlo dpp vs npp data
fit2a = fitlm(table(log(npp0_scat), log(dpp0_scat)), 'linear',...
    'Weights', sqrt(npp_scat)); % fit a linear regression weighted by sqrt of...
        % ...number of primaries

D_2a = fit2a.Coefficients.Estimate(2); % ensemble scaling exponent
k_2a = exp(fit2a.Coefficients.Estimate(1)); % ensemble scaling prefactor

% 95% confidence intervals for scaling properties
ci_2a = coefCI(fit2a);
ci_D_2a = ci_2a(2,:);
ci_k_2a = exp(ci_2a(1,:));

% 95% ci error bars
dci_D_2a = max(ci_D_2a) - D_2a;
dcip_k_2a = max(ci_k_2a) - k_2a;
dcin_k_2a = k_2a - min(ci_k_2a);

% generate the fit data
dpp_fit2a = k_2a * (npp_bc.^D_2a);
ci_dpp_fit2a = [ci_k_2a(1) * (npp_bc.^ci_D_2a(1)),...
    ci_k_2a(2) * (npp_bc.^ci_D_2a(2))];

nexttile(1)

% plot the main fit and CI bounds
plt2a_fit = plot(npp_bc, dpp_fit2a, 'Color', hex2rgb('#374259'),...
    'LineStyle', '-', 'LineWidth', 1.5);
plt2a_err = plot(npp_bc, ci_dpp_fit2a(:,1), npp_bc, ci_dpp_fit2a(:,2),...
    'Color', hex2rgb('#374259'), 'LineStyle', ':', 'LineWidth', 1);

lgd2a_fit = strcat(string(newline), 'Linear regression fit (weighted by $n_\mathrm{pp}$)',...
    string(newline), '$D_\mathrm{fit}$ =', {' '}, num2str(D_2a, '%.3f'),...
    {' '}, '$\pm$', {' '}, num2str(dci_D_2a, '%.3f'), {','}, {' '},...
    '$k_\mathrm{fit}$ =', {' '}, num2str(k_2a, '%.2f'),...
    {' '}, {'$\pm$'}, {' '}, num2str(max(dcip_k_2a, dcin_k_2a), '%.2f'));

legend(cat(2, plt2a_bc, plt2a_scat, plt2a_fit),...
    cat(2, lgd2a_bc, lgd2a_scat, lgd2a_fit),...
    'interpreter', 'latex', 'FontSize', 11, 'location', 'southoutside');

%%% find best fitting to the dpp vs da data
fit2b = fitlm(table(log(da_scat), log(dpp0_scat)), 'linear', 'Weights',...
    sqrt(npp_scat)); % fit a linear regression weighted by sqrt of...
        % ...number of primaries

D_2b = fit2b.Coefficients.Estimate(2); % ensemble scaling exponent
k_2b = exp(fit2b.Coefficients.Estimate(1)); % ensemble scaling prefactor

% 95% confidence intervals for scaling properties
ci_2b = coefCI(fit2b);
ci_D_2b = ci_2b(2,:);
ci_k_2b = exp(ci_2b(1,:));

% 95% ci error bars
dci_D_2b = max(ci_D_2b) - D_2b;
dcip_k_2b = max(ci_k_2b) - k_2b;
dcin_k_2b = k_2b - min(ci_k_2b);

% generate the fit data
dpp_fit2b = k_2b * (da_uc.^D_2b);
ci_dpp_fit2b = [ci_k_2b(1) * (da_uc.^ci_D_2b(1)),...
    ci_k_2b(2) * (da_uc.^ci_D_2b(2))];

nexttile(2)

% plot the main fit and CI bounds
plt2b_fit = plot(da_uc, dpp_fit2b, 'Color', hex2rgb('#632626'),...
    'LineStyle', '-', 'LineWidth', 1.5);
plt2b_err = plot(da_uc, ci_dpp_fit2b(:,1), da_uc, ci_dpp_fit2b(:,2),...
    'Color', hex2rgb('#632626'), 'LineStyle', ':', 'LineWidth', 1);

lgd2b_fit = strcat(string(newline), 'Linear regression fit (weighted by $n_\mathrm{pp}$)', ...
    string(newline), '$D_\mathrm{fit}$ =', {' '}, num2str(D_2b, '%.2f'),...
    {' '}, '$\pm$', {' '}, num2str(dci_D_2b, '%.2f'), {','}, {' '},...
    '$k_\mathrm{fit}$ =', {' '}, num2str(k_2b, '%.2f'),...
    {' '}, {'$\pm$'}, {' '}, num2str(max(dcip_k_2b, dcin_k_2b), '%.2f'));

legend(cat(2, plt2b_uc, plt2b_scat, plt2b_fit),...
    cat(2, lgd2b_uc, lgd2b_scat, lgd2b_fit),...
    'interpreter', 'latex', 'FontSize', 11, 'location', 'southoutside');

%% Correct random seeds, regenerate distribution plots and find best fits %%

% initialize figure
f3 = figure(3);
f3.Position = [150, 150, 900, 600];
set(f3, 'color', 'white');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact')

nexttile(1) % dpp vs. npp figure

% plot combination of Brasil's and universal correlations
plt3a_bc = copyobj(plt2a_bc, f3.CurrentAxes);
hold on

% plot random gussian points with npp rounded to integer numbers
plt3a_scat = scatter(npp_scat, dpp_scat, 12, hex2rgb('#789DBC'), '.',...
    'LineWidth', 1);
lgd3a_scat = strcat(string(newline), 'Gaussian noise seeds with $n_\mathrm{pp}$ rounded');

% find best curve fit
fit3a = fitlm(table(log(npp_scat), log(dpp_scat)), 'linear',...
    'Weights', sqrt(npp_scat)); % fit a linear regression weighted by sqrt of...
        % ...number of primaries

D_3a = fit3a.Coefficients.Estimate(2); % ensemble scaling exponent
k_3a = exp(fit3a.Coefficients.Estimate(1)); % ensemble scaling prefactor

% 95% confidence intervals for scaling properties
ci_3a = coefCI(fit3a);
ci_D_3a = ci_3a(2,:);
ci_k_3a = exp(ci_3a(1,:));

% 95% ci error bars
dci_D_3a = max(ci_D_3a) - D_3a;
dcip_k_3a = max(ci_k_3a) - k_3a;
dcin_k_3a = k_3a - min(ci_k_3a);

% generate the fit data
dpp_fit3a = k_3a * (npp_bc.^D_3a);
ci_dpp_fit3a = [ci_k_3a(1) * (npp_bc.^ci_D_3a(1)),...
    ci_k_3a(2) * (npp_bc.^ci_D_3a(2))];

% plot the main fit and CI bounds
plt3a_fit = plot(npp_bc, dpp_fit3a, 'Color', hex2rgb('#374259'),...
    'LineStyle', '-', 'LineWidth', 1.5);
plt3a_err = plot(npp_bc, ci_dpp_fit3a(:,1), npp_bc, ci_dpp_fit3a(:,2),...
    'Color', hex2rgb('#374259'), 'LineStyle', ':', 'LineWidth', 1);

lgd3a_fit = strcat(string(newline), 'Linear regression fit (weighted by $n_\mathrm{pp}$)',...
    string(newline), '$D_\mathrm{fit}$ =', {' '}, num2str(D_3a, '%.3f'),...
    {' '}, '$\pm$', {' '}, num2str(dci_D_3a, '%.3f'), {','}, {' '},...
    '$k_\mathrm{fit}$ =', {' '}, num2str(k_3a, '%.2f'),...
    {' '}, {'$\pm$'}, {' '}, num2str(max(dcip_k_3a, dcin_k_3a), '%.2f'));

legend(cat(2, plt3a_bc, plt3a_scat, plt3a_fit),...
    cat(2, lgd2a_bc, lgd3a_scat, lgd3a_fit),...
    'interpreter', 'latex', 'FontSize', 11, 'location', 'southoutside');

% appearance configs for dpp vs npp subplot
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([0.8 * min(npp_scat), 1.2 * max(npp_scat)])
ylim([0.9 * min(dpp_scat), 1.1 * max(dpp_scat)])
xlabel('$n_\mathrm{pp}$ [-]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
box on

nexttile(2) % dpp vs. da figure

% plot universal correlation
plt3b_uc = copyobj(plt2b_uc, f3.CurrentAxes);
hold on

% plot library in dpp vs. da domain
plt3b_scat = scatter(da_scat, dpp_scat, 12, hex2rgb('#B06161'), '.',...
    'LineWidth', 1);
lgd3b_scat = strcat(string(newline), 'Transformation of Gaussian noise seeds', string(newline),...
    'from $d_\mathrm{pp}$-$n_\mathrm{pp}$ space to $d_\mathrm{pp}$-$d_\mathrm{a}$ space',...
    string(newline), 'using correlation of Brasil et al. (1999)');

% find best curve fit
fit3b = fitlm(table(log(da_scat), log(dpp_scat)), 'linear',...
    'Weights', sqrt(npp_scat)); % fit a linear regression weighted by sqrt of...
        % ...number of primaries

D_3b = fit3b.Coefficients.Estimate(2); % ensemble scaling exponent
k_3b = exp(fit3b.Coefficients.Estimate(1)); % ensemble scaling prefactor

% 95% confidence intervals for scaling properties
ci_3b = coefCI(fit3b);
ci_D_3b = ci_3b(2,:);
ci_k_3b = exp(ci_3b(1,:));

% 95% ci error bars
dci_D_3b = max(ci_D_3b) - D_3b;
dcip_k_3b = max(ci_k_3b) - k_3b;
dcin_k_3b = k_3b - min(ci_k_3b);

% generate the fit data
dpp_fit3b = k_3b * (da_uc.^D_3b);
ci_dpp_fit3b = [ci_k_3b(1) * (da_uc.^ci_D_3b(1)),...
    ci_k_3b(2) * (da_uc.^ci_D_3b(2))];

% plot the main fit and CI bounds
plt3b_fit = plot(da_uc, dpp_fit3b, 'Color', hex2rgb('#632626'),...
    'LineStyle', '-', 'LineWidth', 1.5);
plt3b_err = plot(da_uc, ci_dpp_fit3b(:,1), da_uc, ci_dpp_fit3b(:,2),...
    'Color', hex2rgb('#632626'), 'LineStyle', ':', 'LineWidth', 1);

lgd3b_fit = strcat(string(newline), 'Linear regression fit (weighted by $n_\mathrm{pp}$)',...
    string(newline), '$D_\mathrm{fit}$ =', {' '}, num2str(D_3b, '%.2f'),...
    {' '}, '$\pm$', {' '}, num2str(dci_D_3b, '%.2f'), {','}, {' '},...
    '$k_\mathrm{fit}$ =', {' '}, num2str(k_3b, '%.2f'),...
    {' '}, {'$\pm$'}, {' '}, num2str(max(dcip_k_3b, dcin_k_3b), '%.2f'));

legend(cat(2, plt3b_uc, plt3b_scat, plt3b_fit),...
    cat(2, lgd2b_uc, lgd3b_scat, lgd3b_fit),...
    'interpreter', 'latex', 'FontSize', 11, 'location', 'southoutside');

% appearance configs for dpp vs da subplot
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([0.9 * min(da_scat), 1.1 * max(da_scat)])
ylim([0.9 * min(dpp_scat), 1.1 * max(dpp_scat)])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
box on

%% select seeds and scale aggregates

% randomize the order of aggregates
ind_rand = randperm(n_agg_raw);
pars_out.pp = cat(1, pars_raw.pp(ind_rand));
pars_out.n = cat(1, pars_raw.n(ind_rand));
pars_out.dpp_g = cat(1, pars_raw.dpp_g(ind_rand,1));
pars_out.da = cat(1, pars_raw.da(ind_rand));

% logical variable for the following pass/reject loop
chk_agg = true(n_agg_raw, 1); % whether an aggregate is previously not...
    % ...selected (i.e. still available)

rpp_scat = zeros(n_scat,1); % initialize scaling ratios

ind_agg_raw = 1 : n_agg_raw; % to keep track of aggregate id selected for scaling

% assign aggregates to scatter seeds 
for i = 1 : n_scat
    
    ii = ind_agg_raw(chk_agg); % first load all aggregate ids available

    % find the closest (to the seed) aggregate in terms of dpp/da ratio 
    ii0 = find(abs((dpp_scat(i) ./ da_scat(i)) -...
        (pars_out.dpp_g(chk_agg) ./ pars_out.da(chk_agg))) ==...
        min(abs((dpp_scat(i) ./ da_scat(i)) -...
        (pars_out.dpp_g(chk_agg) ./ pars_out.da(chk_agg)))), 1);
    
    ii = ii(ii0); % find actual id of aggregate selected
        
    % scale aggregate to the selected seed
    rpp_scat(i) = (1e-9) * dpp_scat(i) / pars_out.dpp_g(ii,1);
    pars_out.pp{ii}(:,2:5) = repmat(rpp_scat(i), pars_out.n(ii), 4) .*...
        pars_out.pp{ii}(:,2:5);
    
    chk_agg(ii) = false; % won't be using this seed another time
    
end

% only keep scaled aggregates
pars_out.pp = cat(1, pars_out.pp(~chk_agg));
pars_out.n = cat(1, pars_out.n(~chk_agg));

% recalculate projected area based on new scaling
pars_out.da = 2 * sqrt(PAR.PROJECTION(pars_out, [], n_mc, n_ang) / pi); % sanity check
% pars_out.da = rpp_scat .* cat(1, pars_out.da(~chk_agg));

% calculate GM and GSD of primary particle size for the scaled aggregates
pars_out = PAR.SIZING(pars_out);

% find projected area size distribution parameters of aggregates...
    % ...after Monte Carlo sampling
GM_da_out = geomean(pars_out.da); % geometric mean
GSD_da_out = UTILS.GEOSTD(pars_out.da); % geometric standard deviation

% find ensemble primary particle size distribution metrics
pp_ens = cat(1, pars_out.pp{:}); % concatinate primaries across aggregates
GM_dpp_ens_out = geomean(pp_ens(:,2)); % geometric mean
GSD_dpp_ens_out = UTILS.GEOSTD(pp_ens(:,2)); % geometric standard deviation

% metrics for primary particle size distribution within the aggregates
GM_dpp_out = geomean(pars_out.dpp_g(:,1)); % geometric mean of...
    % ...geomtric mean of primary particle size within aggregates
GSD_dpp_out = UTILS.GEOSTD(pars_out.dpp_g(:,1)); % geometric standard...
    % ...deviation of geomtric mean of primary particle size within aggregates
mu_sigmapp_out = mean(pars_out.dpp_g(:,2)); % arithmetic mean of...
    % ...geomtric standard deviation of primary particle size within aggregates
sd_sigmapp_out = std(pars_out.dpp_g(:,2)); % arithmetic standard deviation of...
    % ...geomtric standard deviation of primary particle size within aggregates

%% Plot scaled aggregates

% initialize figure
f4= figure(4);
f4.Position = [200, 200, 900, 600];
set(f4, 'color', 'white');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact')

nexttile(1) % dpp vs. npp figure

% plot combination of Brasil's and universal correlations
plt4a_bc = copyobj(plt2a_bc, f4.CurrentAxes);
hold on

% plot scaled aggregates
plt4a_scat = scatter(pars_out.n, 1e9 * pars_out.dpp_g(:,1), 5,...
    hex2rgb('#789DBC'), 'o', 'LineWidth', 1);
lgd4a_scat = strcat(string(newline), 'Aggregates scaled to a random Gaussian',...
    string(newline), 'bivariate distribution');

% find best curve fit
fit4a = fitlm(table(log(pars_out.n), log(1e9 * pars_out.dpp_g(:,1))),...
    'linear', 'Weights', sqrt(pars_out.n)); % fit a linear regression...
        % ...weighted by sqrt of number of primaries

D_4a = fit4a.Coefficients.Estimate(2); % ensemble scaling exponent
k_4a = exp(fit4a.Coefficients.Estimate(1)); % ensemble scaling prefactor

% 95% confidence intervals for scaling properties
ci_4a = coefCI(fit4a);
ci_D_4a = ci_4a(2,:);
ci_k_4a = exp(ci_4a(1,:));

% 95% ci error bars
dci_D_4a = max(ci_D_4a) - D_4a;
dcip_k_4a = max(ci_k_4a) - k_4a;
dcin_k_4a = k_4a - min(ci_k_4a);

% generate the fit data
dpp_fit4a = k_4a * (npp_bc.^D_4a);
ci_dpp_fit4a = [ci_k_4a(1) * (npp_bc.^ci_D_4a(1)),...
    ci_k_4a(2) * (npp_bc.^ci_D_4a(2))];

% plot the main fit and CI bounds
plt4a_fit = plot(npp_bc, dpp_fit4a, 'Color', hex2rgb('#374259'),...
    'LineStyle', '-', 'LineWidth', 1.5);
plt4a_err = plot(npp_bc, ci_dpp_fit4a(:,1), npp_bc, ci_dpp_fit4a(:,2),...
    'Color', hex2rgb('#374259'), 'LineStyle', ':', 'LineWidth', 1);

lgd4a_fit = strcat(string(newline), 'Linear regression fit (weighted by $n_\mathrm{pp}$)',...
    string(newline), '$D_\mathrm{fit}$ =', {' '}, num2str(D_4a, '%.3f'),...
    {' '}, '$\pm$', {' '}, num2str(dci_D_4a, '%.3f'), {','}, {' '},...
    '$k_\mathrm{fit}$ =', {' '}, num2str(k_4a, '%.2f'),...
    {' '}, {'$\pm$'}, {' '}, num2str(max(dcip_k_4a, dcin_k_4a), '%.2f'));

legend(cat(2, plt4a_bc, plt4a_scat, plt4a_fit),...
    cat(2, lgd2a_bc, lgd4a_scat, lgd4a_fit),...
    'interpreter', 'latex', 'FontSize', 11, 'location', 'southoutside');

% appearance configs for dpp vs npp subplot
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([0.8 * min(pars_out.n), 1.2 * max(pars_out.n)])
ylim([1e9 * 0.9 * min(pars_out.dpp_g(:,1)), 1e9 * 1.1 * max(pars_out.dpp_g(:,1))])
xlabel('$n_\mathrm{pp}$ [-]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
box on

nexttile(2) % dpp vs. da figure

% plot universal correlation
plt4b_uc = copyobj(plt2b_uc, f4.CurrentAxes);
hold on

% plot library in dpp vs. da domain
plt4b_scat = scatter(1e9 * pars_out.da, 1e9 * pars_out.dpp_g(:,1), 8,...
    hex2rgb('#B06161'), 'o', 'LineWidth', 1);
lgd4b_scat = strcat(string(newline), 'Aggregates scaled to a random Gaussian',...
    string(newline), 'bivariate distribution');

% find best curve fit
fit4b = fitlm(table(log(1e9 * pars_out.da), log(1e9 * pars_out.dpp_g(:,1))),...
    'linear', 'Weights', sqrt(pars_out.n)); % fit a linear regression...
        % ...weighted by sqrt of number of primaries

D_4b = fit4b.Coefficients.Estimate(2); % ensemble scaling exponent
k_4b = exp(fit4b.Coefficients.Estimate(1)); % ensemble scaling prefactor

% 95% confidence intervals for scaling properties
ci_4b = coefCI(fit4b);
ci_D_4b = ci_4b(2,:);
ci_k_4b = exp(ci_4b(1,:));

% 95% ci error bars
dci_D_4b = max(ci_D_4b) - D_4b;
dcip_k_4b = max(ci_k_4b) - k_4b;
dcin_k_4b = k_3b - min(ci_k_4b);

% generate the fit data
dpp_fit4b = k_4b * (da_uc.^D_4b);
ci_dpp_fit4b = [ci_k_4b(1) * (da_uc.^ci_D_4b(1)),...
    ci_k_4b(2) * (da_uc.^ci_D_4b(2))];

% plot the main fit and CI bounds
plt4b_fit = plot(da_uc, dpp_fit4b, 'Color', hex2rgb('#632626'),...
    'LineStyle', '-', 'LineWidth', 1.5);
plt4b_err = plot(da_uc, ci_dpp_fit4b(:,1), da_uc, ci_dpp_fit4b(:,2),...
    'Color', hex2rgb('#632626'), 'LineStyle', ':', 'LineWidth', 1);

lgd4b_fit = strcat(string(newline), 'Linear regression fit (weighted by $n_\mathrm{pp}$)',...
    string(newline), '$D_\mathrm{fit}$ =', {' '}, num2str(D_4b, '%.2f'),...
    {' '}, '$\pm$', {' '}, num2str(dci_D_4b, '%.2f'), {','}, {' '},...
    '$k_\mathrm{fit}$ =', {' '}, num2str(k_4b, '%.2f'),...
    {' '}, {'$\pm$'}, {' '}, num2str(max(dcip_k_4b, dcin_k_4b), '%.2f'));

legend(cat(2, plt4b_uc, plt4b_scat, plt4b_fit),...
    cat(2, lgd2b_uc, lgd4b_scat, lgd4b_fit),...
    'interpreter', 'latex', 'FontSize', 11, 'location', 'southoutside');

% appearance configs for dpp vs da subplot
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([1e9 * 0.9 * min(pars_out.da(:,1)), 1e9 * 1.1 * max(pars_out.da(:,1))])
ylim([1e9 * 0.9 * min(pars_out.dpp_g(:,1)), 1e9 * 1.1 * max(pars_out.dpp_g(:,1))])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
box on

%% Draw resulting aggregates %%

f5 = figure(5);
f5.Position = [250, 100, 900, 900];
set(f5, 'color', 'white');

tiledlayout(3, 3, 'Padding', 'none', 'TileSpacing', 'none')

jj = sort(randperm(n_scat, 9));

for j = 1 : 9
    nexttile
    UTILS.PLOTPP(pars_out.pp{jj(j)}(:,3), pars_out.pp{jj(j)}(:,4),...
        pars_out.pp{jj(j)}(:,5), pars_out.pp{jj(j)}(:,2))
end

%% Export plots %%

exportgraphics(f1, 'outputs\raw.png',...
    'BackgroundColor','none', 'ContentType','vector', 'Resolution', 300)
exportgraphics(f2, 'outputs\seeds_initial.png',...
    'BackgroundColor','none', 'ContentType','vector', 'Resolution', 300)
exportgraphics(f3, 'outputs\seeds_corrected.png',...
    'BackgroundColor','none', 'ContentType','vector', 'Resolution', 300)
exportgraphics(f4, 'outputs\scaled.png',...
    'BackgroundColor','none', 'ContentType','vector', 'Resolution', 300)
exportgraphics(f5, 'outputs\render.png',...
    'BackgroundColor','none', 'ContentType','vector', 'Resolution', 300)
