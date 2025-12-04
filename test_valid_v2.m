clear; clc; close all;

% load scaled uniform aggregates - Brasil's values
fdir_dpp_da_1 = 'D:\Hamed\CND\PhD\Weekly\2025\Week_of_27OCT25';
fname_dpp_da_1 = 'BiggerRepo-Scat135-Real1-Frac0055-27OCT25';
vars_dpp_da = {'pars_out'};
load(fullfile(fdir_dpp_da_1, fname_dpp_da_1), vars_dpp_da{:});
pars_out_dalpha_108_k_alpha_110 = pars_out;

% load scaled uniform aggregates - D_alpha=1.20
fdir_dpp_da_2 = 'D:\Hamed\CND\PhD\Weekly\2025\Week_of_27OCT25\Results\d_alpha_1.2';
fname_dpp_da_2 = 'BiggerRepo-Scat135-Real1-Frac0055-10NOV25-Dalpha120';
vars_dpp_da = {'pars_out'};
load(fullfile(fdir_dpp_da_2, fname_dpp_da_2), vars_dpp_da{:});
pars_out_dalpha_120_k_alpha_110 = pars_out;

% load scaled uniform aggregates - k_alpha=1.20
fdir_dpp_da_3 = 'D:\Hamed\CND\PhD\Weekly\2025\Week_of_27OCT25\Results\k_alpha_1.2';
fname_dpp_da_3 = 'BiggerRepo-Scat135-Real1-Frac0055-10NOV25-Kalpha120';
vars_dpp_da = {'pars_out'};
load(fullfile(fdir_dpp_da_3, fname_dpp_da_3), vars_dpp_da{:});
pars_out_dalpha_108_k_alpha_120 = pars_out;

clear pars_out

% second-stage langevin dynamics data to be imported
fdir_LD2 = 'D:\Hamed\CND\PhD\Weekly\2025\Week_of_27OCT25\Results\Brasil\LD2-28-Oct-2025_BiggerRepo-Scat135-Real1-Frac0055-27OCT25';
fname_LD2 = 'LD2_30-Oct-2025_12-50-54_Final';
varnames_LD2 = {'fl'};
load(fullfile(fdir_LD2, fname_LD2), varnames_LD2{:});

% load experimental data for validation
fdir_exp = 'D:\Hamed\CND\PhD\Publication\Experiment\Effective-Density-Compiled';
fname_exp = 'Effective-Density-Compiled_17-Mar-2025_05-00-44';
load(strcat(fdir_exp, '\', fname_exp, '.mat'), 'dist_grp')

% calculate mobility diameter for all aggregates
dm{1} = TRANSP.DIAMOBIL(pars_out_dalpha_108_k_alpha_110.dg,...
    pars_out_dalpha_108_k_alpha_110.da, fl);
dm{2} = TRANSP.DIAMOBIL(pars_out_dalpha_120_k_alpha_110.dg,...
    pars_out_dalpha_120_k_alpha_110.da, fl);
dm{3} = TRANSP.DIAMOBIL(pars_out_dalpha_108_k_alpha_120.dg,...
    pars_out_dalpha_108_k_alpha_120.da, fl);

% calculate their effective density
rho_eff{1} = cellfun(@(x) 1860 * (pi/6) * sum(x(:,2).^3),...
    pars_out_dalpha_108_k_alpha_110.pp) ./ dm{1}.^3;
rho_eff{2} = cellfun(@(x) 1860 * (pi/6) * sum(x(:,2).^3),...
    pars_out_dalpha_120_k_alpha_110.pp) ./ dm{2}.^3;
rho_eff{3} = cellfun(@(x) 1860 * (pi/6) * sum(x(:,2).^3),...
    pars_out_dalpha_108_k_alpha_120.pp) ./ dm{3}.^3;

% convert mobility diameter to [nm]
dm{1} = 1e9*dm{1};
dm{2} = 1e9*dm{2};
dm{3} = 1e9*dm{3};

% initialize figure
f = figure(1);
f.Position = [100, 100, 450, 500];
set(f, 'color', 'white');

% draw universal correlation (for effective density)
n_uc = 100; % number of data points to draw the relation
D_m = 2.48; % exponent
rho_eff_100 = 510; % pefactor
dm_lim_uc = [min(cat(1,dm{1},dm{2},dm{3})) max(cat(1,dm{1},dm{2},dm{3}))];
r_uc = (dm_lim_uc(2) / dm_lim_uc(1)) ^ (1 / (n_uc - 1));
dm_uc = dm_lim_uc(1) * ones(n_uc,1) .* r_uc .^ (((1 : n_uc) - 1)');
rho_eff_uc = rho_eff_100 * (dm_uc/100) .^ (D_m-3);
plt_uc = plot(dm_uc, rho_eff_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

% appearance configs for plot
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 16)
ylabel('$\rho_\mathrm{eff}$ [nm]', 'interpreter', 'latex', 'FontSize', 16)
box on

%%% comparing bayesian fits for effective density trends

% initialize strcuture for bayesian fits
bayesfit = struct('xfit', cell(3,1), 'yfit', cell(3,1),...
    'bounds_yfit', cell(3,1), 'afit', cell(3,1), 'bounds_afit', cell(3,1));

% colors for fits
colors = [0.7 0.9 0.7; 0.9 0.7 0.7; 0.7 0.7 0.9]; % pastel green, red, blue

% initialize plot placeholder
plt = cell(3,1);

for i = 1:3
    
    % derive the fit
    [bayesfit(i).yfit, bayesfit(i).xfit, bayesfit(i).bounds_yfit,...
        bayesfit(i).afit, bayesfit(i).bounds_afit] =...
        UTILS.BAYESFIT_POLY2(dm{i}, rho_eff{i}, ones(length(dm{i}),1), 1000);
    
    % plot the fit
    plt{i} = loglog(bayesfit(i).xfit, bayesfit(i).yfit, 'Color', colors(i,:),...
        'LineWidth', 2);

    % confidence interval
    fill([bayesfit(i).xfit; flipud(bayesfit(1).xfit)],...
        [bayesfit(i).bounds_yfit(:,1); flipud(bayesfit(i).bounds_yfit(:,2))],...
        colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);

end

plt_valid = scatter(dist_grp(1).d_mode, dist_grp(1).rho_eff,...
    20, hex2rgb('#C96868'), 'v', 'LineWidth', 1.5);

% generate legend
legend([plt_uc; cat(1,plt{:}); plt_valid], {'Olfert $\&$ Rogak (2019)',...
    'D=1.08 k=1.10', 'D=1.20 k=1.10', 'D=1.08 k=1.20', 'Experiment'},...
    'interpreter', 'latex', 'FontSize', 12, 'Location', 'northoutside',...
    'NumColumns', 2, 'Orientation', 'horizontal');


