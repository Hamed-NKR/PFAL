clear; clc; close all;

%% loading data %%

% import simulation results from second-stage langevin dynamics
% fdir_LD2 = 'D:\Hamed\CND\PhD\Publication\Paper2\Library_Final\1_35\LD2';
% fname_LD2 = 'LD2-temp';
fdir_LD2 = 'D:\Hamed\CND\PhD\Publication\Paper2\Validation\Scatter\LD2-25-Sep-2025_Pre_LD2_Scatter___25_09_25';
fname_LD2 = 'LD2_25-Sep-2025_13-03-12_Final';
varnames_LD2 = {'parsdata', 'fl'};
load(fullfile(fdir_LD2, fname_LD2), varnames_LD2{:});
ind_sim = [1 4];
rows_to_merge = 4;

% import experimental effective density data from tandem measurements
fdir_exp = 'D:\Hamed\CND\PhD\Publication\Experiment\Effective-Density-Compiled';
fname_exp = 'Effective-Density-Compiled_17-Mar-2025_05-00-44';
load(strcat(fdir_exp, '\', fname_exp, '.mat'), 'dist_grp')
ind_exp = [1 3];

%% preprocess data

% decide which rows in parsdata to merge with the row after them
parsdata = merge_pars_rows(parsdata, rows_to_merge);

% only keep simulations/experiments of interest
parsdata = parsdata(ind_sim); 
dist_grp = dist_grp(ind_exp);

%% 2. Validating effective density %%

% number of simulation/experimental datasets
n_sim = length(ind_sim); n_exp = length(ind_exp);

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

% initialize figure
f2 = figure(2);
f2.Position = [150, 150, 450, 500];
set(f2, 'color', 'white');

% draw universal correlation (for effective density)
n_uc = 100; % number of data points to draw the relation
D_m = 2.48; % exponent
rho_eff_100 = 510; % pefactor
dm_lim_uc = [0.9 * min(min(cat(1,dm_sim{:})), min(cat(2, dist_grp(:).d_mode))),...
    1.2 * max(max(cat(1,dm_sim{:})), max(cat(2, dist_grp(:).d_mode)))];
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
bayesfit = struct('xfit', cell(n_sim,1), 'yfit', cell(n_sim,1),...
    'bounds_yfit', cell(n_sim,1), 'afit', cell(n_sim,1),...
    'bounds_afit', cell(n_sim,1));

% colors for fits
colors_sim = [hex2rgb('#F9C0AF'); hex2rgb('#C4DFE5')];
colors_exp = [hex2rgb('#C66E52'); hex2rgb('#758A93')];

% initialize plot & legend placeholders
plt_sim = cell(n_sim,1); lgd_sim = cell(n_sim,1);
plt_exp = cell(n_exp,1); lgd_exp = cell(n_exp,1);

for i = 1 : min(length(ind_sim), length(ind_exp))
    
    % derive the bayesian fit on the popul;ation of simulated aggregates
    [bayesfit(i).yfit, bayesfit(i).xfit, bayesfit(i).bounds_yfit,...
        bayesfit(i).afit, bayesfit(i).bounds_afit] =...
        UTILS.BAYESFIT_POLY2(dm_sim{i}, rho_eff_sim{i}, ones(length(dm_sim{i}),1), 1000);
    
    % plot the fit
    plt_sim{i} = loglog(bayesfit(i).xfit, bayesfit(i).yfit, 'Color', colors_sim(i,:),...
        'LineWidth', 2);

    lgd_sim{i} = sprintf('Simulation %d', i); % update legend

    % confidence interval
    fill([bayesfit(i).xfit; flipud(bayesfit(i).xfit)],...
        [bayesfit(i).bounds_yfit(:,1); flipud(bayesfit(i).bounds_yfit(:,2))],...
        colors_sim(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    % draw experimental datapoints
    plt_exp{i} = scatter(dist_grp(i).d_mode, dist_grp(i).rho_eff,...
    20, colors_exp(i,:), 'v', 'LineWidth', 1.5);
    
    lgd_exp{i} = sprintf('Experiment %d', i); % update legend

end

% generate legend
legend([plt_uc; cat(1,plt_sim{:}); cat(1,plt_exp{:})],...
    [{'Olfert $\&$ Rogak (2019)'}; lgd_sim; lgd_exp],...
    'interpreter', 'latex', 'FontSize', 12, 'Location', 'northoutside',...
    'NumColumns', 2, 'Orientation', 'horizontal');

% adjust axes
xlim([0.9 * max(min(cat(1,dm_sim{:})), min(cat(2, dist_grp(:).d_mode))),...
    1.1 * max(max(cat(1,dm_sim{:})), max(cat(2, dist_grp(:).d_mode)))])
ylim([0.9 * min(min(cat(1,rho_eff_sim{:})), min(cat(2, dist_grp(:).rho_eff))),...
    1.1 * max(max(cat(1,rho_eff_sim{:})), max(cat(2, dist_grp(:).rho_eff)))])

%%% Possible causes for deviation in effective density %%%
%   ** uncertainty in material density (i.e. 1860 kg/m3)
%   ** bias in calibration coefficients of AAC
%   ** neglecting overlapping of primary particles
%   ** uncertainty in manual sizing of primary particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

function parsdata = merge_pars_rows(parsdata, rows_to_merge)
% "MERGE_PARS_ROWS" merges selected rows of parsdata with the row after
%   them and remove duplicate aggregates.
%
%   - parsdata      : struct array with fields dpp, da, n_hyb, pp, etc.
%   - rows_to_merge : vector of indices (in the CURRENT parsdata) such that
%                     each i in rows_to_merge is merged with i+1.
%
%   Example: rows_to_merge = [1 3] merges row 1 with 2 and row 3 with 4.

    % work from the bottom up so index shifting doesn't break things
    rows_to_merge = sort(rows_to_merge(:), 'descend');

    for kk = 1:numel(rows_to_merge)
        i = rows_to_merge(kk);
        j = i + 1;

        if j > numel(parsdata)
            error('merge_pars_rows:badIndex', ...
                  'Row %d has no row after it to merge with.', i);
        end

        % concatenate aggregates from rows i and j
        merged = parsdata(i); % start from row i

        % assuming these fields exist in parsdata; add more if needed
        merged.dpp   = [parsdata(i).dpp; parsdata(j).dpp];
        merged.sigmapp   = [parsdata(i).sigmapp; parsdata(j).sigmapp];
        merged.da    = [parsdata(i).da; parsdata(j).da];
        merged.dm   = [parsdata(i).dm; parsdata(j).dm];
        merged.dg   = [parsdata(i).dg; parsdata(j).dg];
        merged.n_hyb = [parsdata(i).n_hyb; parsdata(j).n_hyb];
        merged.pp    = [parsdata(i).pp; parsdata(j).pp];
        merged.npp    = [parsdata(i).npp; parsdata(j).npp];

        nagg = numel(merged.pp);

        % find duplicates based on primary-particle IDs
        ind_flt = false(nagg,1);

        for a = 1:nagg-1
            if ind_flt(a), continue; end
            pa = sort(unique(merged.pp{a}(:,1)));

            for b = a+1:nagg
                if ind_flt(b), continue; end
                pb = sort(unique(merged.pp{b}(:,1)));

                if isequal(pa, pb)
                    ind_flt(b) = true; % mark as duplicate
                end
            end
        end

        keep = ~ind_flt;

        % apply filter to each per-aggregate field
        merged.dpp = merged.dpp(keep);
        merged.sigmapp = merged.sigmapp(keep);
        merged.da = merged.da(keep);
        merged.dm = merged.dm(keep);
        merged.dg = merged.dg(keep);
        merged.n_hyb = merged.n_hyb(keep);
        merged.pp = merged.pp(keep);
        merged.npp = merged.npp(keep);

        % put merged row back into parsdata
        parsdata(i) = merged; % replace row i with merged version
        parsdata(j) = []; % delete row j

    end
end
