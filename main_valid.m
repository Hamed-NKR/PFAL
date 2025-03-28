clc
clear
clf('reset')
close all

%% initialization

% address of simulation data to be imported
fdir_simul_0 = 'C:\Users\hmdnkr\Documents\GitHub\MCEM\outputs\postLD2-27-Mar-2025_17-58-33_LD2-25NOV24';
fname_simul_0 = 'Post_LD2-25NOV24';

% variables of interest in the simulation data
varnames_simul = {'parsdata', 'dm', 'rho_eff', 'pars_ens'};

ii0 = [1,3,4]; % data group indices in simulation data to be plotted

% address of experimental data to be improted (for validation)
fdir_exp = 'F:\Experiment\Effective-Density-Compiled';
fname_exp = 'Effective-Density-Compiled_17-Mar-2025_05-00-44';

%% figure for temporal distribution of simulation effective density vs. experiments

% load the simulation to be validated
for i = 1 : numel(varnames_simul)
    load(strcat(fdir_simul_0, '\', fname_simul_0, '.mat'),...
        varnames_simul{i})
end

% rename imported variables
parsdata0 = parsdata;
dm0 = dm;
rho_eff_0 = rho_eff;
pars_ens_0 = pars_ens;

clear parsdata dm rho_eff pars_ens % clear original names

% initialize figure
f1 = figure(2);
f1.Position = [50, 50, 900, 500];
set(f1, 'color', 'white')

% initialize layout
tl1 = tiledlayout(1,2);
tl1.TileSpacing = 'compact';
tl1.Padding = 'compact';

% initialize placholders for plots & legends
plt11 = cell(6, 1);
legtxt11 = cell(6, 1);
plt12 = cell(5, 1);
legtxt12 = cell(5, 1);

% initialize marker colors, shapes and sizes
mc1 = colormap(hot);
cind2 = round(1 + (length(mc1) - 1) .* (0.05 : 0.7 / (5 - 1) : 0.75)');
mc1 = mc1(cind2,:);
ms1 = [8, 16, 12, 16, 8];
mt1 = {'^', 's', 'p', '*', 'o'};
mc1(end,:) = [236,230,61] / 255;

nexttile(1)

% universal correlation for effective density vs. mobility diameter
D_m = 2.48; % exponent
rho_eff_100 = 510; % pefactor
dm_lim_uc = [1e0 2e4];  % limits on the mobility diameter
n_dm_uc = 1e4; % number of data
uc2 = @(y) rho_eff_100 * (y / 100) .^ (D_m - 3); % on-demand function...
    % ...for the universal correlation
r_uc2 = (dm_lim_uc(2) / dm_lim_uc(1)) ^ (1 / (n_dm_uc - 1)); % increments
dm_uc = dm_lim_uc(1) * ones(n_dm_uc,1) .* r_uc2 .^ (((1 : n_dm_uc) - 1)'); % mobility setpoints
rho_eff_uc = uc2(dm_uc); % universal effective densities
plt11{end} = plot(dm_uc, rho_eff_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on
legtxt11{end} = 'Olfert $\&$ Rogak (2019)';

rho0 = 1860; % material density for soot

for i = ii0
    
    ii = find(i==ii0,1);

    plt11{ii} = scatter(1e9 * dm0{i}, rho_eff_0{i},...
        ms1(i), mc1(i,:), mt1{i}, 'LineWidth', 1);
    
    if i == 1
        legtxt11(ii) = strcat('$n_\mathrm{agg}/n_\mathrm{agg_0}$ =',...
            {' '}, num2str(parsdata(i).r_n_agg(1), '%.0f'));
    elseif i == 3
        legtxt11(ii) = strcat('$n_\mathrm{agg}/n_\mathrm{agg_0}$ =',...
            {' '}, num2str(parsdata(i).r_n_agg(1), '%.1f'));
    else
        legtxt11(ii) = strcat('$n_\mathrm{agg}/n_\mathrm{agg_0}$ =',...
            {' '}, num2str(parsdata(i).r_n_agg(1), '%.2f'));
    end
end

bounds_dm_f1 = [1e9 * 0.9 * min(cat(1, dm0{ii0})),...
    1e9 * 1.1 * max(cat(1, dm0{ii0}))];
bounds_rho_f1 = [0.9 * min(cat(1, rho_eff_0{ii0})),...
    1.1 * max(cat(1, rho_eff_0{ii0}))];

% load experimental data
load(strcat(fdir_exp, '\', fname_exp, '.mat'), 'dist_grp')

% plot experimental data
plt11{4} = scatter(dist_grp(1).d_mode, dist_grp(1).rho_eff,...
    25, hex2rgb('#C96868'), 'v', 'LineWidth', 1);
legtxt11{4} = 'Lo-Aglom';
plt11{5} = scatter(dist_grp(3).d_mode, dist_grp(3).rho_eff,...
    35, hex2rgb('#7EACB5'), 'h', 'LineWidth', 1);
legtxt11{5} = 'Hi-Aglom';

% set plot appearances
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim(bounds_dm_f1)
ylim(bounds_rho_f1)
xlabel('$d_\mathrm{m} [nm]$', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$\rho_\mathrm{eff} [kg/m^3]$', 'interpreter', 'latex',...
    'FontSize', 14)
legend(cat(1, plt11{:}), legtxt11, 'interpreter', 'latex',...
    'FontSize', 11, 'NumColumns', 2, 'Location', 'southoutside')
