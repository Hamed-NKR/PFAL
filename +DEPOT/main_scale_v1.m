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

% metrics of projected area size distribution
gm_da = 83.22; % geometric mean
gsd_da = 1.67; % geometric standard deviation

gm_dpp = 16.67; % ensemble geometric mean of primary particle diameter distribution
gsd_dpp = 1.22; % ensemble geometric standard deviation

cn_scale = 0.6; % proportion of random aggregates chosen for filtering

% address of aggregate library to be imported for scaling and filtering
fdir = 'D:\Hamed\CND\PhD\My Articles\DLCA2\mainscatter_sigmapp13';
fname = 'wsp-sigmapp-13';
varname = 'pp0';
vardir = '';

% resolution for projected area calculation
n_mc = 1e4;
n_ang = 25;

n_bin_filter = 20;

%% Raw data against the Brasil's and universal correlations %%

% load non-scaled first-stage library data
load(strcat(fdir, '\', fname, '.mat'), varname)

% create particle structure  
pars_raw.pp = eval(strcat(varname, vardir)); % store primary particle info

% correct the structure if necessary
if logical(nnz(size(pars_raw.pp)))
    
    % Convert the pp data into a 1d cell array
    pars_raw.pp = pars_raw.pp(:);
    
    % Remove unused cells
    pars_raw.pp = pars_raw.pp(~cellfun('isempty', pars_raw.pp));
    
    % Merge pp info from different times
    pars_raw.pp = cat(1, pars_raw.pp{:});

end

n_agg_raw = length(pars_raw.pp); % initial number of aggregates
pars_raw.n = zeros(n_agg_raw,1);
for i = 1 : n_agg_raw
    % count the number of primaries for each aggregate
    pars_raw.n(i) = size(pars_raw.pp{i}, 1); 
end
pars_raw = PAR.SIZING(pars_raw); % get characteristic sizes

eval(['clear ', varname]) % delete pp0

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

%% Scale and filter aaggregates to the universal correlation %%

% number of random aggregates to be selected after filtering
n_scale = round(cn_scale * length(pars_raw.n));

% make a structure for scaling the aggregates
pars_scale.n = cat(1, pars_raw.n);
pars_scale.pp = cat(1, pars_raw.pp);

% get target mean primary particle sizes for scaling 
dpp_scale = buc(pars_scale.n);

% get the scaling factors
rpp_scale = (1e-9) * dpp_scale ./ cat(1, pars_raw.dpp_g(:,1));

% implement scaling on each aggregate
for i = 1 : length(pars_scale.n)
    pars_scale.pp{i}(:,2:5) = repmat(rpp_scale(i), pars_scale.n(i), 4) .*...
        pars_scale.pp{i}(:,2:5) * gm_dpp / uc(gm_da);
end

% recalculate projected area based on new scaling
da_scale = 2 * sqrt(PAR.PROJECTION(pars_scale, [], n_mc, n_ang) / pi);        

% filter aggregates through a lognormal distribution on the projected area
[da_filter, ind_filter] = TRANSP.LNSAMPLING(da_scale, 1e-9 * gm_da,...
    gsd_da, n_scale, n_bin_filter);

% Remove empty size bins
da_filter = da_filter(~cellfun('isempty', da_filter)); 
ind_filter = ind_filter(~cellfun('isempty', ind_filter));

% Compile data from different size bins
da_filter = cat(1, da_filter{:}); 
ind_filter = cat(1, ind_filter{:});

% store the filtered aggregates and make the output aggregate structure
pars_out.pp = pars_scale.pp(ind_filter);
pars_out.n = pars_scale.n(ind_filter);
pars_out.da = da_filter;

% calculate GM and GSD of primary particle size for the scaled aggregates
pars_out = PAR.SIZING(pars_out);

% find projected area size distribution parameters of aggregates...
    % ...after Monte Carlo sampling
GM_da_out = geomean(cat(1, pars_out.da)); % geometric mean
GSD_da_out = UTILS.GEOSTD(cat(1, pars_out.da)); % geometric standard deviation

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

% initialize figure for output scaled+filtered aggregates
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

% plot scaled+filtered aggregates in dpp vs npp domain
plt2a_scale = scatter(pars_out.n, 1e9 * pars_out.dpp_g(:,1), 15,...
    hex2rgb('#789DBC'), 'o', 'LineWidth', 1);
lgd2a_scale = strcat(string(newline), 'Transformation of simulated aggregates',...
    string(newline), 'to follow correlations of Brasil et al. (1999)',...
    string(newline), 'and Olfert \& Rogak (2019)');

% appearance configs for dpp vs npp subplot
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([0.8 * min(pars_out.n), 1.2 * max(pars_out.n)])
ylim([0.9 * 1e9 * min(pars_out.dpp_g(:,1)),...
    1.1 * 1e9 * max(pars_out.dpp_g(:,1))])
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
plt2b_scale = scatter(1e9 * pars_out.da, 1e9 * pars_out.dpp_g(:,1), 15,...
    hex2rgb('#B06161'), 'o', 'LineWidth', 1);
lgd2b_scale = strcat(string(newline), 'Transformation of simulated aggregates',...
    string(newline), 'to follow correlations of Brasil et al. (1999)',...
    string(newline), 'and Olfert \& Rogak (2019)');

% appearance configs for dpp vs da subplot
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([0.9 * 1e9 * min(pars_out.da), 1.1 * 1e9 * max(pars_out.da)])
ylim([0.9 * 1e9 * min(pars_out.dpp_g(:,1)),...
    1.1 * 1e9 * max(pars_out.dpp_g(:,1))])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
box on

%% evaluate trends of "non-rounded" Monte Carlo points  %%

figure(f2);

%%% find the best curve fit to the Monte Carlo dpp vs npp data
fit2a = fitlm(table(log(pars_out.n), log(1e9 * pars_out.dpp_g(:,1))),...
    'linear', 'Weights', sqrt(pars_out.n)); % fit a linear regression...
        % ...weighted by sqrt of number of primaries

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

legend(cat(2, plt2a_bc, plt2a_scale, plt2a_fit),...
    cat(2, lgd2a_bc, lgd2a_scale, lgd2a_fit),...
    'interpreter', 'latex', 'FontSize', 11, 'location', 'southoutside');

%%% find best fitting to the dpp vs da data
fit2b = fitlm(table(log(1e9 * pars_out.da), log(1e9 * pars_out.dpp(:,1))),...
    'linear', 'Weights', sqrt(pars_out.n)); % fit a linear regression...
        % ...weighted by sqrt of number of primaries

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

legend(cat(2, plt2b_uc, plt2b_scale, plt2b_fit),...
    cat(2, lgd2b_uc, lgd2b_scale, lgd2b_fit),...
    'interpreter', 'latex', 'FontSize', 11, 'location', 'southoutside');

%% Draw resulting aggregates %%

f3 = figure(3);
f3.Position = [250, 100, 900, 900];
set(f3, 'color', 'white');

tiledlayout(3, 3, 'Padding', 'none', 'TileSpacing', 'none')

jj = sort(randperm(n_scale, 9));

for j = 1 : 9
    nexttile
    UTILS.PLOTPP(pars_out.pp{jj(j)}(:,3), pars_out.pp{jj(j)}(:,4),...
        pars_out.pp{jj(j)}(:,5), pars_out.pp{jj(j)}(:,2))
end

%% Export plots %%

exportgraphics(f1, 'outputs\raw_distribution.png',...
    'BackgroundColor','none', 'ContentType','vector', 'Resolution', 300)
exportgraphics(f2, 'outputs\scaled_filtered_distribution.png',...
    'BackgroundColor','none', 'ContentType','vector', 'Resolution', 300)
exportgraphics(f3, 'outputs\renders.png',...
    'BackgroundColor','none', 'ContentType','vector', 'Resolution', 300)
