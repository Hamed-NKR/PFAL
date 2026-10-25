clear; clc; close all;

%% loading data %%

% import simulation results from second-stage langevin dynamics
fdir_LD2 = 'D:\Hamed\CND\PhD\Publication\Paper2\Library_Final\1_35\LD2\Continuum';
fname_LD2 = 'LD2_11-Dec-2025_19-13-29_Final';
varnames_LD2 = {'parsdata', 'fl'};
load(fullfile(fdir_LD2, fname_LD2), varnames_LD2{:});
ind_sim = [1 4];
rows_to_merge = 4;

% import experimental effective density data from tandem measurements
fdir_exp = 'D:\Hamed\CND\PhD\Publication\Experiment\Effective-Density-Compiled';
fname_exp = 'Effective-Density-Compiled_17-Mar-2025_05-00-44';
load(strcat(fdir_exp, '\', fname_exp, '.mat'), 'dist_grp')
ind_exp = [1 3];

% import manually processed TEM image data
fdir_tem = 'D:\Hamed\CND\PhD\Publication\Paper2\TEM_Manual';
fname_tem = 'TEM_Manual';
varnames_tem = {'Aggs_lal_1', 'Aggs_exdil',...
    'dbarpp_manu_lal_1', 'dbarpp_manu_exdil', 'id_agg_exdil'};
load(fullfile(fdir_tem, fname_tem), varnames_tem{:});

% number of simulation/experimental datasets
n_sim = length(ind_sim); n_exp = length(ind_exp);

% organize TEM data
da_tem = {cat(1, Aggs_lal_1.da); cat(1, Aggs_exdil.da)};
dpp_tem = {dbarpp_manu_lal_1; dbarpp_manu_exdil};
ipp = {1:length(da_tem{1}); id_agg_exdil}; % indices for aggregates whose...
    % ...primary particles have been manually sized
clear 'Aggs_lal_1' 'Aggs_exdil' 'dbarpp_manu_lal_1' 'dbarpp_manu_exdil'...
    'id_agg_exdil' % for tidiness


%% preprocessing data %%

% decide which rows in parsdata to merge with the row after them
parsdata = UTILS.MERGE_PARSDATA_ROWS(parsdata, rows_to_merge);

% only keep simulations/experiments of interest
parsdata = parsdata(ind_sim); 
dist_grp = dist_grp(ind_exp);

%% initializing graphical appearance %%

% set up appearance configs for fits and datapoints
colors_sim_fit = [hex2rgb('#B77466'); hex2rgb('#758A93')];
colors_sim_agg = [hex2rgb('#E16A54'); hex2rgb('#434E78')];
colors_exp = [hex2rgb('#8D493A'); hex2rgb('#606676')];
markerSize_sim = [4 , 8];
markerSymbol_exp = {'v' , 'o'};
markerSize_exp = [15 , 15];
lgd_txt = {'Lo-Aglom', 'Hi-Aglom'};

%% 1. Validating primary particle diameter vs. projected area diameter %%

% initialize figure 1
f1 = figure(1);
f1.Position = [50, 50, 450, 500];
set(f1, 'color', 'white');

% draw universal correlation (in dpp-da domain)
n1_uc = 100; % number of data points to draw the relation
D_TEM = 0.35; % exponent
dpp_100 = 17.8; % pefactor
da_lim_uc = [0.9 * min(1e9*min(cat(1,parsdata.da)), min(cat(1, da_tem{:}))),...
    1.2 * max(1e9*max(cat(1,parsdata.da)), max(cat(1, da_tem{:})))];
r1_uc = (da_lim_uc(2) / da_lim_uc(1)) ^ (1 / (n1_uc - 1));
da_uc = da_lim_uc(1) * ones(n1_uc,1) .* r1_uc .^ (((1 : n1_uc) - 1)');
dpp_uc = dpp_100 * (da_uc/100) .^ D_TEM;
plt1_uc = plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

% set up axes for plot
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 16)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 16)
box on

% initialize strcuture for bayesian fits to dpp-da data cloud
bayesfit1 = struct('xfit', cell(n_sim,1), 'yfit', cell(n_sim,1),...
    'bounds_yfit', cell(n_sim,1), 'afit', cell(n_sim,1),...
    'bounds_afit', cell(n_sim,1));

% initialize plot & legend placeholders
plt1_sim = cell(n_sim,1); lgd_sim = cell(n_sim,1);
plt1_exp = cell(n_exp,1); lgd_exp = cell(n_exp,1);

for i = 1 : min(length(ind_sim), length(ind_exp))
    
    % draw raw simulated data first (background)
    scatter(1e9*parsdata(i).da, 1e9*parsdata(i).dpp, markerSize_sim(i),...
        colors_sim_agg(i,:), 'filled',...
        'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);

    % derive the bayesian fit on the population of simulated aggregates
    [bayesfit1(i).yfit, bayesfit1(i).xfit, bayesfit1(i).bounds_yfit,...
        bayesfit1(i).afit, bayesfit1(i).bounds_afit] =...
            UTILS.FIT_POLY(1e9*parsdata(i).da, 1e9*parsdata(i).dpp,...
            ones(length(parsdata(i).da),1), 1000,'Degree', 1,...
            'XTransform', 'log10', 'YTransform', 'log10',...
            'SlopeOffset', 0);
    
    % plot the fit
    plt1_sim{i} = loglog(bayesfit1(i).xfit, bayesfit1(i).yfit,...
        'Color', colors_sim_fit(i,:), 'LineWidth', 2);

    lgd_sim{i} = sprintf('Simulation: %s', lgd_txt{i}); % update legend

    % confidence interval
    fill([bayesfit1(i).xfit; flipud(bayesfit1(i).xfit)],...
        [bayesfit1(i).bounds_yfit(:,1); flipud(bayesfit1(i).bounds_yfit(:,2))],...
        colors_sim_fit(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    % draw experimental datapoints
    plt1_exp{i} = scatter(da_tem{i}(ipp{i}), dpp_tem{i},...
    markerSize_exp(i), colors_exp(i,:), markerSymbol_exp{i},...
    'LineWidth', 1.5);
    
    lgd_exp{i} = sprintf('Experiment: %s', lgd_txt{i}); % update legend

end

% generate legend
legend([cat(1,plt1_sim{:}); cat(1,plt1_exp{:}); plt1_uc],...
    [lgd_sim; lgd_exp; {'Olfert $\&$ Rogak (2019)'}],...
    'interpreter', 'latex', 'FontSize', 11, 'Location', 'northoutside',...
    'NumColumns', 2, 'Orientation', 'horizontal');

% adjust axes
xlim([0.95 * min(1e9*min(cat(1,parsdata.da)), min(cat(1,da_tem{:}))),...
    1.05 * max(1e9*max(cat(1,parsdata.da)),...
    max(cat(1,da_tem{1}(ipp{1}),da_tem{2}(ipp{2}))))])
ylim([0.95 * min(1e9*min(cat(1,parsdata.dpp)), min(cat(1,dpp_tem{:}))),...
    1.05 * max(1e9*max(cat(1,parsdata.dpp)), max(cat(1,dpp_tem{:})))])


%% 2. Validating effective density vs. mobility diameter %%

% initialize mobility diameter and effective density for simulated...
    % ...aggregates
dm_sim = cell(n_sim,1);
rho_eff_sim = cell(n_sim,1);

% calculate mobility diameter for simulated aggregates
for i = 1 : length(ind_sim)
    
    % pars{i}.pp = parsdata(i).pp;
    % pars{i}.n = parsdata(i).npp;
    % da{i} = 2 * sqrt(PAR.PROJECTION(pars{i}, [], 1e3, 10) / pi);
    % dm{i} = TRANSP.DIAMOBIL(parsdata(i).dg, da{i}, fl);
    dm_sim{i} = TRANSP.DIAMOBIL(parsdata(i).dg, parsdata(i).da, fl);
    
    % calculate their effective density
    rho_eff_sim{i} = cellfun(@(x) 1860 * (pi/6) * sum(x(:,2).^3),...
        parsdata(i).pp) ./ dm_sim{i}.^3;
    
    % convert mobility diameter to [nm]
    dm_sim{i} = 1e9*dm_sim{i};

end

% initialize figure 2
f2 = figure(2);
f2.Position = [150, 150, 450, 500];
set(f2, 'color', 'white');

% draw universal correlation (for effective density)
n2_uc = 100; % number of data points to draw the relation
D_m = 2.48; % exponent
rho_eff_100 = 510; % pefactor
dm_lim_uc = [0.9 * min(min(cat(1,dm_sim{:})), min(cat(2, dist_grp(:).d_mode))),...
    1.2 * max(max(cat(1,dm_sim{:})), max(cat(2, dist_grp(:).d_mode)))];
r2_uc = (dm_lim_uc(2) / dm_lim_uc(1)) ^ (1 / (n2_uc - 1));
dm_uc = dm_lim_uc(1) * ones(n2_uc,1) .* r2_uc .^ (((1 : n2_uc) - 1)');
rho_eff_uc = rho_eff_100 * (dm_uc/100) .^ (D_m-3);
plt2_uc = plot(dm_uc, rho_eff_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

% set up axes for plot
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 16)
ylabel('$\rho_\mathrm{eff}$ [kg/$m^3$]', 'interpreter', 'latex', 'FontSize', 16)
box on

% initialize strcuture for bayesian fits for effective density trends
bayesfit2 = struct('xfit', cell(n_sim,1), 'yfit', cell(n_sim,1),...
    'bounds_yfit', cell(n_sim,1), 'afit', cell(n_sim,1),...
    'bounds_afit', cell(n_sim,1));

% initialize plot & legend placeholders
plt2_sim = cell(n_sim,1);
plt2_exp = cell(n_exp,1);

for i = 1 : min(length(ind_sim), length(ind_exp))
    
    % draw raw simulated data first (background)
    scatter(dm_sim{i}, rho_eff_sim{i}, markerSize_sim(i),...
        colors_sim_agg(i,:), 'filled',...
        'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);

    % derive the bayesian fit on the population of simulated aggregates
    [bayesfit2(i).yfit, bayesfit2(i).xfit, bayesfit2(i).bounds_yfit,...
        bayesfit2(i).afit, bayesfit2(i).bounds_afit] =...
            UTILS.FIT_POLY(dm_sim{i}, rho_eff_sim{i},...
            ones(length(dm_sim{i}),1), 1000,'Degree', 2,...
            'XTransform', 'log10', 'YTransform', 'log10',...
            'SlopeOffset', 3);
    
    % plot the fit
    plt2_sim{i} = loglog(bayesfit2(i).xfit, bayesfit2(i).yfit,...
        'Color', colors_sim_fit(i,:), 'LineWidth', 2);
    
    % confidence interval
    fill([bayesfit2(i).xfit; flipud(bayesfit2(i).xfit)],...
        [bayesfit2(i).bounds_yfit(:,1); flipud(bayesfit2(i).bounds_yfit(:,2))],...
        colors_sim_fit(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    % draw experimental datapoints
    plt2_exp{i} = scatter(dist_grp(i).d_mode, dist_grp(i).rho_eff,...
    markerSize_exp(i), colors_exp(i,:), markerSymbol_exp{i},...
    'LineWidth', 1.5);
    
end

% generate legend
legend([cat(1,plt2_sim{:}); cat(1,plt2_exp{:}); plt2_uc],...
    [lgd_sim; lgd_exp; {'Olfert $\&$ Rogak (2019)'}],...
    'interpreter', 'latex', 'FontSize', 11, 'Location', 'northoutside',...
    'NumColumns', 2, 'Orientation', 'horizontal');

% adjust axes
xlim([0.95 * min(min(cat(1,dm_sim{:})), min(cat(2, dist_grp(:).d_mode))),...
    1.05 * max(max(cat(1,dm_sim{:})), max(cat(2, dist_grp(:).d_mode)))])
ylim([0.95 * min(min(cat(1,rho_eff_sim{:})), min(cat(2, dist_grp(:).rho_eff))),...
    1.05 * max(max(cat(1,rho_eff_sim{:})), max(cat(2, dist_grp(:).rho_eff)))])

% ----- Possible causes for deviation in effective density ----- %
%   (1) uncertainty in material density (i.e. 1860 kg/m3)
%   (2) bias in calibration coefficients of AAC
%   (3) neglecting overlapping of primary particles
%   (4) uncertainty in manual sizing of primary particles
%   (5) selecting modes from the tandem SMPS size distributions
% -------------------------------------------------------------- %

