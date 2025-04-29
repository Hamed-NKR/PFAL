clc
clear
clf('reset')
close all

%% initialization %%

% address of simulation data to be imported
fdir_simul_13_scat = 'F:\DLCA2\outputs\postLD2-28-Mar-2025_18-40-55_LD2-25NOV24';
fname_simul_13_scat = 'Post_LD2-25NOV24';
fdir_simul_13_flat = 'F:\DLCA2\outputs\postLD2-28-Mar-2025_18-36-43_LD2_27-Nov-2024_04-21-58_Final';
fname_simul_13_flat = 'Post_LD2_27-Nov-2024_04-21-58_Final';
fdir_simul_10_scat = 'F:\DLCA2\outputs\postLD2-28-Mar-2025_18-57-18_LD2-25-Nov-2024_19-37-48_Final';
fname_simul_10_scat = 'Post_LD2-25-Nov-2024_19-37-48_Final';
fdir_simul_10_flat = 'F:\DLCA2\outputs\postLD2-28-Mar-2025_18-59-48_LD2-26-Nov-2024_23-56-35_Final';
fname_simul_10_flat = 'Post_LD2-26-Nov-2024_23-56-35_Final';

% variables of interest in the simulation data
varnames_simul = {'parsdata', 'pars_flt', 'bayesfit'};

ii0 = [1,3,4]; % data group indices in simulation data to be plotted

% address of experimental data to be improted (for validation)
fdir_exp = 'F:\Experiment\Effective-Density-Compiled';
fname_exp = 'Effective-Density-Compiled_17-Mar-2025_05-00-44';

bayesresol = 500; % number of increment for bayesian fit 

%% load all simulations %%

% load the simulation (gammapp = 1.3, scattered correlation)
for i = 1 : numel(varnames_simul)
    load(strcat(fdir_simul_13_scat, '\', fname_simul_13_scat, '.mat'),...
        varnames_simul{i})
end
% rename imported variables
parsdata_13_scat = parsdata; pars_flt_13_scat = pars_flt;
bayesfit_13_scat = bayesfit;
clear parsdata pars_flt bayesfit % clear original names

% gammapp = 1.3, flat correlation
for i = 1 : numel(varnames_simul)
    load(strcat(fdir_simul_13_flat, '\', fname_simul_13_flat, '.mat'),...
        varnames_simul{i})
end
parsdata_13_flat = parsdata; pars_flt_13_flat = pars_flt;
bayesfit_13_flat = bayesfit;
clear parsdata pars_flt bayesfit

% gammapp = 1.0, scattered correlation
for i = 1 : numel(varnames_simul)
    load(strcat(fdir_simul_10_scat, '\', fname_simul_10_scat, '.mat'),...
        varnames_simul{i})
end
parsdata_10_scat = parsdata; pars_flt_10_scat = pars_flt;
bayesfit_10_scat = bayesfit;
clear parsdata pars_flt bayesfit

% gammapp = 1.0, flat correlation
for i = 1 : numel(varnames_simul)
    load(strcat(fdir_simul_10_flat, '\', fname_simul_10_flat, '.mat'),...
        varnames_simul{i})
end
parsdata_10_flat = parsdata; pars_flt_10_flat = pars_flt;
bayesfit_10_flat = bayesfit;
clear parsdata pars_flt bayesfit

%% temporal distribution of simulation effective density vs. experiments %%

% initialize figure
f1 = figure(1);
f1.Position = [50, 50, 1100, 650];
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

for i = ii0
    
    ii = find(i==ii0,1);

    plt11{ii} = scatter(1e9 * parsdata_13_scat(i).dm, parsdata_13_scat(i).rho_eff,...
        ms1(i), mc1(i,:), mt1{i}, 'LineWidth', 1);
    
    if i == 1
        legtxt11(ii) = strcat('$n_\mathrm{agg}/(n_\mathrm{agg})_2$ =',...
            {' '}, num2str(parsdata_13_scat(i).r_n_agg(1), '%.0f'));
    elseif i == 3
        legtxt11(ii) = strcat('$n_\mathrm{agg}/(n_\mathrm{agg})_2$ =',...
            {' '}, num2str(parsdata_13_scat(i).r_n_agg(1), '%.1f'));
    else
        legtxt11(ii) = strcat('$n_\mathrm{agg}/(n_\mathrm{agg})_2$ =',...
            {' '}, num2str(parsdata_13_scat(i).r_n_agg(1), '%.2f'));
    end
end

bounds_dm_f1 = [1e9 * 0.9 * min(cat(1, parsdata_13_scat(ii0).dm)),...
    1e9 * 1.1 * max(cat(1, parsdata_13_scat(ii0).dm))];
bounds_rho_f1 = [0.9 * min(cat(1, parsdata_13_scat(ii0).rho_eff)),...
    1.1 * max(cat(1, parsdata_13_scat(ii0).rho_eff))];

% load experimental data
load(strcat(fdir_exp, '\', fname_exp, '.mat'), 'dist_grp')

% plot experimental data
plt11{4} = scatter(dist_grp(1).d_mode, dist_grp(1).rho_eff,...
    30, hex2rgb('#C96868'), 'v', 'LineWidth', 1.5);
legtxt11{4} = 'Lo-Aglom';
plt11{5} = scatter(dist_grp(3).d_mode, dist_grp(3).rho_eff,...
    45, hex2rgb('#7EACB5'), 'h', 'LineWidth', 1.5);
legtxt11{5} = 'Hi-Aglom';

% set plot appearances
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim(bounds_dm_f1)
ylim(bounds_rho_f1)
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 18)
ylabel('$\rho_\mathrm{eff} \mathrm{[kg/m^3]}$', 'interpreter', 'latex',...
    'FontSize', 18)
legend(cat(1, plt11{:}), legtxt11, 'interpreter', 'latex',...
    'FontSize', 14, 'NumColumns', 2, 'Location', 'southoutside')

nexttile(2)

% replot universal correlation for second tile
plt12{end} = plot(dm_uc, rho_eff_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on
legtxt12{end} = 'Olfert $\&$ Rogak (2019)';

% initialize strcuture for bayesian fits to effective density validation data
bayesfit_valid = struct('xfit', cell(4,1), 'yfit', cell(4,1),...
    'bounds_yfit', cell(4,1), 'afit', cell(4,1), 'bounds_afit', cell(4,1));

% bayesian fit to non-agglomerated simulated aggregates
[bayesfit_valid(1).yfit, bayesfit_valid(1).xfit, bayesfit_valid(1).bounds_yfit,...
    bayesfit_valid(1).afit, bayesfit_valid(1).bounds_afit] =...
    UTILS.BAYESFIT_POLY4(1e9 * parsdata_13_scat(1).dm,...
    parsdata_13_scat(1).rho_eff, parsdata_13_scat(1).r_n_agg, bayesresol);
plt12{1} = loglog(bayesfit_valid(1).xfit, bayesfit_valid(1).yfit, 'Color', mc1(1,:),...
    'LineWidth', 2);
fill([bayesfit_valid(1).xfit; flipud(bayesfit_valid(1).xfit)],...
    [bayesfit_valid(1).bounds_yfit(:,1); flipud(bayesfit_valid(1).bounds_yfit(:,2))],...
    mc1(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
legtxt12{1} = legtxt11{1};
bayesfit_valid(1).name = legtxt12{1};

% bayesian fit to highly-agglomerated simulated aggregates...
    % ...(0.03 < r_n_agg < 0.1)
[bayesfit_valid(2).yfit, bayesfit_valid(2).xfit, bayesfit_valid(2).bounds_yfit,...
    bayesfit_valid(2).afit, bayesfit_valid(2).bounds_afit] =...
    UTILS.BAYESFIT_POLY4(1e9 * cat(1, parsdata_13_scat(3:4).dm),...
    cat(1, parsdata_13_scat(3:4).rho_eff), cat(1, parsdata_13_scat(3:4).r_n_agg),...
    bayesresol);
plt12{2} = loglog(bayesfit_valid(2).xfit, bayesfit_valid(2).yfit, 'Color', mc1(4,:),...
    'LineWidth', 2);
fill([bayesfit_valid(2).xfit; flipud(bayesfit_valid(2).xfit)],...
    [bayesfit_valid(2).bounds_yfit(:,1); flipud(bayesfit_valid(2).bounds_yfit(:,2))],...
    mc1(4,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
legtxt12(2) = strcat(num2str(parsdata_13_scat(3).r_n_agg(1), '%.1f'), {' '},...
    '$\leq n_\mathrm{agg}/(n_\mathrm{agg})_2 \leq$', {' '},...
    num2str(parsdata_13_scat(4).r_n_agg(1), '%.2f'));
bayesfit_valid(2).name = legtxt12{2};

% bayesian fit to 'Lo-Aglom' experimental condition
[bayesfit_valid(3).yfit, bayesfit_valid(3).xfit, bayesfit_valid(3).bounds_yfit,...
    bayesfit_valid(3).afit, bayesfit_valid(3).bounds_afit] =...
    UTILS.BAYESFIT_POLY4(dist_grp(1).d_mode', dist_grp(1).rho_eff',...
    ones(size(dist_grp(1).rho_eff')), bayesresol);
plt12{3} = loglog(bayesfit_valid(3).xfit, bayesfit_valid(3).yfit, 'Color',...
    hex2rgb('#C96868'), 'LineWidth', 2.5, 'LineStyle', ':');
fill([bayesfit_valid(3).xfit; flipud(bayesfit_valid(3).xfit)],...
    [bayesfit_valid(3).bounds_yfit(:,1); flipud(bayesfit_valid(3).bounds_yfit(:,2))],...
    hex2rgb('#C96868'), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
legtxt12{3} = 'Lo-Aglom';
bayesfit_valid(3).name = legtxt12{3};

% bayesian fit to 'Hi-Aglom' experimental condition
[bayesfit_valid(4).yfit, bayesfit_valid(4).xfit, bayesfit_valid(4).bounds_yfit,...
    bayesfit_valid(4).afit, bayesfit_valid(4).bounds_afit] =...
    UTILS.BAYESFIT_POLY4(dist_grp(3).d_mode', dist_grp(3).rho_eff',...
    ones(size(dist_grp(3).rho_eff')), bayesresol);
plt12{4} = loglog(bayesfit_valid(4).xfit, bayesfit_valid(4).yfit, 'Color',...
    hex2rgb('#7EACB5'), 'LineWidth', 2.5, 'LineStyle', ':');
fill([bayesfit_valid(4).xfit; flipud(bayesfit_valid(4).xfit)],...
    [bayesfit_valid(4).bounds_yfit(:,1); flipud(bayesfit_valid(4).bounds_yfit(:,2))],...
    hex2rgb('#7EACB5'), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
legtxt12{4} = 'Hi-Aglom';
bayesfit_valid(4).name = legtxt12{4};

% set plot appearances
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim(bounds_dm_f1)
ylim(bounds_rho_f1)
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 18)
ylabel('$\rho_\mathrm{eff} \mathrm{[kg/m^3]}$', 'interpreter', 'latex',...
    'FontSize', 18)

legend(cat(1, plt12{:}), legtxt12, 'interpreter', 'latex',...
    'FontSize', 14, 'NumColumns', 2, 'Location', 'southoutside')

%% parametric studies on simulations (temporal effective densities) %%

% initialize figure
f2 = figure(2);
f2.Position = [100, 100, 850, 900];
set(f2, 'color', 'white')

% initialize layout
tl2 = tiledlayout(2,2);
tl2.TileSpacing = 'compact';
tl2.Padding = 'compact';

% initialize placholders for plots & legends
plt2 = cell(6,1);
legtxt2 = cell(6,1);

ms2 = [10, 20, 15, 20, 10]; % set marker size

% calculate bounds for subplots
bounds_dm_f2 = [1e9 * 0.9 * min([min(cat(1, parsdata_13_scat.dm)),...
    min(cat(1, parsdata_13_flat.dm)), min(cat(1, parsdata_10_scat.dm)),...
    min(cat(1, parsdata_10_flat.dm))]),...
    1e9 * 1.1 * max([max(cat(1, parsdata_13_scat.dm)),...
    max(cat(1, parsdata_13_flat.dm)), max(cat(1, parsdata_10_scat.dm)),...
    max(cat(1, parsdata_10_flat.dm))])];
bounds_rho_f2 = [0.9 * min([min(cat(1, parsdata_13_scat.rho_eff)),...
    min(cat(1, parsdata_13_flat.rho_eff)), min(cat(1, parsdata_10_scat.rho_eff)),...
    min(cat(1, parsdata_10_flat.rho_eff))]),...
    1.1 * max([max(cat(1, parsdata_13_scat.rho_eff)),...
    max(cat(1, parsdata_13_flat.rho_eff)), max(cat(1, parsdata_10_scat.rho_eff)),...
    max(cat(1, parsdata_10_flat.rho_eff))])];

% make titles for subplots
titxt2 = {'$(\gamma_\mathrm{pp})_1$ = 1.3 -- bivariate',...
    '$(\gamma_\mathrm{pp})_1$ = 1.3 -- univariate',...
    '$(\gamma_\mathrm{pp})_1$ = 1.0 -- bivariate',...
    '$(\gamma_\mathrm{pp})_1$ = 1.0 -- univariate'};

for i = 1 : 4

    nexttile(i)

    % plot universal correlation of of Olfert & Rogak (2019)
    if i == 1
        plt2{6} = plot(dm_uc, rho_eff_uc, 'Color', [0.4940 0.1840 0.5560],...
            'LineStyle', '-.', 'LineWidth', 3);
        legtxt2{6} = 'Olfert $\&$ Rogak (2019)';
    else
        plot(dm_uc, rho_eff_uc, 'Color', [0.4940 0.1840 0.5560],...
            'LineStyle', '-.', 'LineWidth', 3);
    end
    hold on

    for j = 1 : 5

        switch i

            % aggregates with initial polydispersity, and with scatter...
                % ...around universal correlation
            case 1
                plt2{j} = scatter(1e9 * parsdata_13_scat(j).dm,...
                    parsdata_13_scat(j).rho_eff, ms2(j), mc1(j,:),...
                    mt1{j}, 'LineWidth', 1);
                if j == 1
                    legtxt2(j) = strcat('$n_\mathrm{agg}/(n_\mathrm{agg})_2$ =',...
                        {' '}, num2str(parsdata_13_scat(j).r_n_agg(1), '%.0f'));
                elseif ismember(j, [2,3])
                    legtxt2(j) = strcat('$n_\mathrm{agg}/(n_\mathrm{agg})_2$ =',...
                        {' '}, num2str(parsdata_13_scat(j).r_n_agg(1), '%.1f'));
                elseif ismember(j, [4,5])
                    legtxt2(j) = strcat('$n_\mathrm{agg}/(n_\mathrm{agg})_2$ =',...
                        {' '}, num2str(parsdata_13_scat(j).r_n_agg(1), '%.2f'));
                end

            % aggregates with initial polydispersity, but without scatter...
                % ...around universal correlation
            case 2
                    scatter(1e9 * parsdata_13_flat(j).dm,...
                        parsdata_13_flat(j).rho_eff, ms2(j), mc1(j,:),...
                        mt1{j}, 'LineWidth', 1);

            % aggregates without initial polydispersity, but with scatter...
                % ...around universal correlation
            case 3
                    scatter(1e9 * parsdata_10_scat(j).dm,...
                        parsdata_10_scat(j).rho_eff, ms2(j), mc1(j,:),...
                        mt1{j}, 'LineWidth', 1);

            % aggregates without initial polydispersity, and without scatter...
                % ...around universal correlation
            case 4
                    scatter(1e9 * parsdata_10_flat(j).dm,...
                        parsdata_10_flat(j).rho_eff, ms2(j), mc1(j,:),...
                        mt1{j}, 'LineWidth', 1);
        end
    end

    % set plot appearances
    box on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
        'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
    xlim(bounds_dm_f2)
    ylim(bounds_rho_f2)
    title(titxt2{i}, 'interpreter', 'latex', 'FontSize', 14)

end

lgd2 = legend(cat(1, plt2{:}), legtxt2, 'interpreter', 'latex',...
    'FontSize', 14, 'NumColumns', 3);
lgd2.Layout.Tile = 'north';
xlabel(tl2, '$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 18)
ylabel(tl2, '$\rho_\mathrm{eff} \mathrm{[kg/m^3]}$', 'interpreter',...
    'latex', 'FontSize', 18)

%% curvefits to parametric studies %%

% initialize figure
f3 = figure(3);
f3.Position = [150, 150, 900, 500];
set(f3, 'color', 'white')

% initialize layout
tl3 = tiledlayout(1,2);
tl3.TileSpacing = 'compact';
tl3.Padding = 'compact';

% initialize placholders for plots & legends
plt3 = cell(6, 2);
legtxt3 = cell(6, 1);

nexttile(1)

% effective density comparsion
plt3{5,1} = plot(dm_uc, rho_eff_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on
legtxt3{5} = 'Olfert $\&$ Rogak (2019)';

%%% plot universal correlation of of Olfert & Rogak (2019)

% aggregates with initial polydispersity and with scatter around...
    % ...universal correlation
plt3{1,1} = plot(bayesfit_13_scat.xfit, bayesfit_13_scat.yfit,...
    'Color', hex2rgb('#C96868'), 'LineWidth', 1.5);
fill([bayesfit_13_scat.xfit; flipud(bayesfit_13_scat.xfit)],...
    [bayesfit_13_scat.bounds_yfit(:,1);...
    flipud(bayesfit_13_scat.bounds_yfit(:,2))],...
    hex2rgb('#C96868'), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
legtxt3{1} = '$(\gamma_\mathrm{pp})_1$ = 1.3, bivariate';

% aggregates with initial polydispersity, but without scatter around...
    % ...universal correlation
plt3{2,1} = plot(bayesfit_13_flat.xfit, bayesfit_13_flat.yfit,...
    'Color', hex2rgb('#FCDC94'), 'LineWidth', 1.5);
fill([bayesfit_13_flat.xfit; flipud(bayesfit_13_flat.xfit)],...
    [bayesfit_13_flat.bounds_yfit(:,1);...
    flipud(bayesfit_13_flat.bounds_yfit(:,2))],...
    hex2rgb('#FCDC94'), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
legtxt3{2} = '$(\gamma_\mathrm{pp})_1$ = 1.3, univariate';

% aggregates without initial polydispersity, but with scatter around...
    % ...universal correlation
plt3{3,1} = plot(bayesfit_10_scat.xfit, bayesfit_10_scat.yfit,...
    'Color', hex2rgb('#7286D3'), 'LineWidth', 1.5);
fill([bayesfit_10_scat.xfit; flipud(bayesfit_10_scat.xfit)],...
    [bayesfit_10_scat.bounds_yfit(:,1);...
    flipud(bayesfit_10_scat.bounds_yfit(:,2))],...
    hex2rgb('#7286D3'), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
legtxt3{3} = '$(\gamma_\mathrm{pp})_1$ = 1.0, bivariate';

% aggregates without initial polydispersity and without scatter around...
    % ...universal correlation
plt3{4,1} = plot(bayesfit_10_flat.xfit, bayesfit_10_flat.yfit,...
    'Color', hex2rgb('#55AD9B'), 'LineWidth', 1.5);
fill([bayesfit_10_flat.xfit; flipud(bayesfit_10_flat.xfit)],...
    [bayesfit_10_flat.bounds_yfit(:,1);...
    flipud(bayesfit_10_flat.bounds_yfit(:,2))],...
    hex2rgb('#55AD9B'), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
legtxt3{4} = '$(\gamma_\mathrm{pp})_1$ = 1.0, univariate';

% set plot appearances
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
bounds_dm_f3 = [min([min(bayesfit_13_scat.xfit),...
    min(bayesfit_13_flat.xfit), min(bayesfit_10_scat.xfit),...
    min(bayesfit_10_flat.xfit)]), max([max(bayesfit_13_scat.xfit),...
    max(bayesfit_13_flat.xfit), max(bayesfit_10_scat.xfit),...
    max(bayesfit_10_flat.xfit)])];
xlim(bounds_dm_f3)
bounds_rho_f3 = [min([min(bayesfit_13_scat.yfit),...
    min(bayesfit_13_flat.yfit), min(bayesfit_10_scat.yfit),...
    min(bayesfit_10_flat.yfit)]), max([max(bayesfit_13_scat.yfit),...
    max(bayesfit_13_flat.yfit), max(bayesfit_10_scat.yfit),...
    max(bayesfit_10_flat.yfit)])];
ylim(bounds_rho_f3)
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 18)
ylabel('$\rho_\mathrm{eff} \mathrm{[kg/m^3]}$', 'interpreter', 'latex',...
    'FontSize', 18)

% mass-mobility exponent comparison
nexttile(2)
plot(bayesfit_13_scat.xfit, bayesfit_13_scat.afit,...
    bayesfit_13_flat.xfit, bayesfit_13_flat.afit,...
    bayesfit_10_scat.xfit, bayesfit_10_scat.afit,...
    bayesfit_10_flat.xfit, bayesfit_10_flat.afit)

% plot exponent of universal correlation
plt3{5,2} = plot(dm_uc, D_m * ones(size(dm_uc)), 'Color',...
    [0.4940 0.1840 0.5560], 'LineStyle', '-.', 'LineWidth', 3);
hold on

% plot exponent limit for diffusion-limited cluster aggregation
plt3{6,2} = plot(dm_uc, 1.78 * ones(size(dm_uc)),...
    'Color', [0, 0, 0], 'LineStyle', ':', 'LineWidth', 3);
legtxt3{6} = 'Fractal dimension';

% aggregates with initial polydispersity and with scatter around...
    % ...universal correlation
plt3{1,2} = plot(bayesfit_13_scat.xfit, bayesfit_13_scat.afit,...
    'Color', hex2rgb('#C96868'), 'LineWidth', 1.5);
fill([bayesfit_13_scat.xfit; flipud(bayesfit_13_scat.xfit)],...
    [bayesfit_13_scat.bounds_afit(:,1);...
    flipud(bayesfit_13_scat.bounds_afit(:,2))],...
    hex2rgb('#C96868'), 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% aggregates with initial polydispersity, but without scatter around...
    % ...universal correlation
plt3{2,2} = plot(bayesfit_13_flat.xfit, bayesfit_13_flat.afit,...
    'Color', hex2rgb('#FCDC94'), 'LineWidth', 1.5);
fill([bayesfit_13_flat.xfit; flipud(bayesfit_13_flat.xfit)],...
    [bayesfit_13_flat.bounds_afit(:,1);...
    flipud(bayesfit_13_flat.bounds_afit(:,2))],...
    hex2rgb('#FCDC94'), 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% aggregates without initial polydispersity, but with scatter around...
    % ...universal correlation
plt3{3,2} = plot(bayesfit_10_scat.xfit, bayesfit_10_scat.afit,...
    'Color', hex2rgb('#7286D3'), 'LineWidth', 1.5);
fill([bayesfit_10_scat.xfit; flipud(bayesfit_10_scat.xfit)],...
    [bayesfit_10_scat.bounds_afit(:,1);...
    flipud(bayesfit_10_scat.bounds_afit(:,2))],...
    hex2rgb('#7286D3'), 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% aggregates without initial polydispersity and without scatter around...
    % ...universal correlation
plt3{4,2} = plot(bayesfit_10_flat.xfit, bayesfit_10_flat.afit,...
    'Color', hex2rgb('#55AD9B'), 'LineWidth', 1.5);
fill([bayesfit_10_flat.xfit; flipud(bayesfit_10_flat.xfit)],...
    [bayesfit_10_flat.bounds_afit(:,1);...
    flipud(bayesfit_10_flat.bounds_afit(:,2))],...
    hex2rgb('#55AD9B'), 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% set plot appearances
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim(bounds_dm_f3)
bounds_Dm_f3 = [min([min(bayesfit_13_scat.afit),...
    min(bayesfit_13_flat.afit), min(bayesfit_10_scat.afit),...
    min(bayesfit_10_flat.afit)]), max([max(bayesfit_13_scat.afit),...
    max(bayesfit_13_flat.afit), max(bayesfit_10_scat.afit),...
    max(bayesfit_10_flat.afit)])];
ylim(bounds_Dm_f3)
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 18)
ylabel('$D_\mathrm{m}$ [-]', 'interpreter', 'latex', 'FontSize', 18)

lgd3 = legend(cat(1, plt3{:,2}), legtxt3, 'interpreter', 'latex',...
    'FontSize', 14, 'NumColumns', 3);
lgd3.Layout.Tile = 'south';

%% temporal primary particle size vs aggregate size %%

% initialize figure
f4 = figure(4);
f4.Position = [200, 200, 850, 900];
set(f4, 'color', 'white')

% initialize layout
tl4 = tiledlayout(2,2);
tl4.TileSpacing = 'compact';
tl4.Padding = 'compact';

plt4 = cell(6,1); % initialize placholders for plots

% calculate bounds for subplots
bounds_da_f4 = [1e9 * 0.9 * min([min(cat(1, parsdata_13_scat.da)),...
    min(cat(1, parsdata_13_flat.da)), min(cat(1, parsdata_10_scat.da)),...
    min(cat(1, parsdata_10_flat.da))]),...
    1e9 * 1.1 * max([max(cat(1, parsdata_13_scat.da)),...
    max(cat(1, parsdata_13_flat.da)), max(cat(1, parsdata_10_scat.da)),...
    max(cat(1, parsdata_10_flat.da))])];
bounds_dpp_f4 = [1e9 * 0.9 * min([min(cat(1, parsdata_13_scat.dpp)),...
    min(cat(1, parsdata_13_flat.dpp)), min(cat(1, parsdata_10_scat.dpp)),...
    min(cat(1, parsdata_10_flat.dpp))]),...
    1e9 * 1.1 * max([max(cat(1, parsdata_13_scat.dpp)),...
    max(cat(1, parsdata_13_flat.dpp)), max(cat(1, parsdata_10_scat.dpp)),...
    max(cat(1, parsdata_10_flat.dpp))])];

% assign variables for universal correlation
D_TEM = 0.35; % exponent
dpp_100 = 17.8; % pefactor
da_lim_uc = [1e0 2e4];  % limits on the projected area diameter
n_da_uc = 1e4; % number of data
uc1 = @(y) dpp_100 * (y / 100) .^ D_TEM; % on-demand function for the...
    % ...forward correlation in the geometrical domain (dpp as a...
    % ...function of da in [nm])
r_uc1 = (da_lim_uc(2) / da_lim_uc(1)) ^ (1 / (n_da_uc - 1));
da_uc = da_lim_uc(1) * ones(n_da_uc,1) .* r_uc1 .^ (((1 : n_da_uc) - 1)');
dpp_uc = uc1(da_uc);

for i = 1 : 4

    nexttile(i)

    % plot universal correlation of of Olfert & Rogak (2019)
    if i == 1
        plt4{6} = plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
            'LineStyle', '-.', 'LineWidth', 3);
        legtxt2{6} = 'Olfert $\&$ Rogak (2019)';
    else
        plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
            'LineStyle', '-.', 'LineWidth', 3);
    end
    hold on

    for j = 1 : 5

        switch i

            % aggregates with initial polydispersity, and with scatter...
                % ...around universal correlation
            case 1
                plt4{j} = scatter(1e9 * parsdata_13_scat(j).da,...
                    1e9 * parsdata_13_scat(j).dpp, ms2(j), mc1(j,:),...
                    mt1{j}, 'LineWidth', 1);

            % aggregates with initial polydispersity, but without scatter...
                % ...around universal correlation
            case 2
                    scatter(1e9 * parsdata_13_flat(j).da,...
                        1e9 * parsdata_13_flat(j).dpp, ms2(j), mc1(j,:),...
                        mt1{j}, 'LineWidth', 1);

            % aggregates without initial polydispersity, but with scatter...
                % ...around universal correlation
            case 3
                    scatter(1e9 * parsdata_10_scat(j).da,...
                        1e9 * parsdata_10_scat(j).dpp, ms2(j), mc1(j,:),...
                        mt1{j}, 'LineWidth', 1);

            % aggregates without initial polydispersity, and without scatter...
                % ...around universal correlation
            case 4
                    scatter(1e9 * parsdata_10_flat(j).da,...
                        1e9 * parsdata_10_flat(j).dpp, ms2(j), mc1(j,:),...
                        mt1{j}, 'LineWidth', 1);
        end
    end

    % set plot appearances
    box on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
        'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
    xlim(bounds_da_f4)
    ylim(bounds_dpp_f4)
    title(titxt2{i}, 'interpreter', 'latex', 'FontSize', 14)

end

lgd4 = legend(cat(1, plt4{:}), legtxt2, 'interpreter', 'latex',...
    'FontSize', 14, 'NumColumns', 3);
lgd4.Layout.Tile = 'north';
xlabel(tl4, '$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 18)
ylabel(tl4, '$d_\mathrm{pp}$ [nm]', 'interpreter',...
    'latex', 'FontSize', 18)

%% a revised version of simulations against experiments (for thesis) %%

% initialize figure
f5 = figure(5);
f5.Position = [250, 50, 600, 700];
set(f5, 'color', 'white')

% initialize placholders for plots & legends
plt5 = cell(5, 1);
legtxt5 = cell(5, 1);

% plot universal correlation for effective density vs. mobility diameter
plt5{end} = plot(dm_uc, rho_eff_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

%%% plot fits to the simulation data %%%
% no agglomeration (r_n_agg = 1; initial moment in second LD stage)
plt5{1} = loglog(bayesfit_valid(1).xfit, bayesfit_valid(1).yfit,...
    'Color', hex2rgb('#9F5255'), 'LineWidth', 2); 
fill([bayesfit_valid(1).xfit; flipud(bayesfit_valid(1).xfit)],...
    [bayesfit_valid(1).bounds_yfit(:,1); flipud(bayesfit_valid(1).bounds_yfit(:,2))],...
    hex2rgb('#9F5255'), 'EdgeColor', 'none', 'FaceAlpha', 0.15);
% high agglomeration (0.03 < r_n_agg < 0.1)
plt5{2} = loglog(bayesfit_valid(2).xfit, bayesfit_valid(2).yfit,...
    'Color', hex2rgb('#006A71'), 'LineWidth', 2);
fill([bayesfit_valid(2).xfit; flipud(bayesfit_valid(2).xfit)],...
    [bayesfit_valid(2).bounds_yfit(:,1); flipud(bayesfit_valid(2).bounds_yfit(:,2))],...
    hex2rgb('#006A71'), 'EdgeColor', 'none', 'FaceAlpha', 0.15);

% plot experimental data
plt5{3} = scatter(dist_grp(1).d_mode, dist_grp(1).rho_eff,...
    30, hex2rgb('#E16A54'), 'v', 'LineWidth', 1.5); % low agglomeration
plt5{4} = scatter(dist_grp(3).d_mode, dist_grp(3).rho_eff,...
    45, hex2rgb('#9ACBD0'), 'h', 'LineWidth', 1.5); % high agglomeration

% make legends
legtxt5 = {'No-Aglom simulation', 'Hi-Aglom simulation',...
    'Lo-Aglom experiment', 'Hi-Aglom experiment',...
    'Olfert $\&$ Rogak (2019)'};

% set plot appearances
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([25 1000])
ylim([35 1250])
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 18)
ylabel('$\rho_\mathrm{eff} \mathrm{[kg/m^3]}$', 'interpreter', 'latex',...
    'FontSize', 18)
legend(cat(1, plt5{:}), legtxt5, 'interpreter', 'latex',...
    'FontSize', 14, 'NumColumns', 2, 'orientation', 'horizontal',...
    'Location', 'southoutside')

%% save plots and workspace %%

% make a directory to save outputs
dir0_out = datestr(datetime('now'));
dir0_out = regexprep(dir0_out, ':', '-');
dir0_out = regexprep(dir0_out, ' ', '_');
dir_out = strcat('outputs\', 'Valid_', dir0_out, '\');
if ~isfolder(dir_out)
    mkdir(dir_out); % if it doesn't exist, create the directory
end

% save worksapce
save(strcat(dir_out, 'Valid_', dir0_out, '.mat'))

% print figures
exportgraphics(f1, strcat(dir_out, 'simul-vs-exp-v1.jpg'),...
    'BackgroundColor','none', 'Resolution', 300)
exportgraphics(f2, strcat(dir_out, 'simul-param-rho-vs-dm.jpg'),...
    'BackgroundColor','none', 'Resolution', 300)
exportgraphics(f3, strcat(dir_out, 'bayesfit-param-rho-and-Dm.jpg'),...
    'BackgroundColor','none', 'Resolution', 300)
exportgraphics(f4, strcat(dir_out, 'simul-param-dpp-vs-da.jpg'),...
    'BackgroundColor','none', 'Resolution', 300)
exportgraphics(f5, strcat(dir_out, 'simul-vs-exp-v2.jpg'),...
    'BackgroundColor','none', 'Resolution', 300)

