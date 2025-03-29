clc
clear
clf('reset')
close all
warning('off')

%% Initialize and load second-stage LD aggregates %%

% address of second-stage langevin dynamics data to be imported
% fdir_simul = 'F:\DLCA2\mainscatter_sigmapp13\SCAT';
% fname_simul = 'LD2-25NOV24';
% fdir_simul = 'F:\DLCA2\mainscatter_sigmapp13\FLAT';
% fname_simul = 'LD2_27-Nov-2024_04-21-58_Final';
% fdir_simul = 'F:\DLCA2\mainscatter_sigmapp10\SCAT';
% fname_simul = 'LD2-25-Nov-2024_19-37-48_Final';
fdir_simul = 'F:\DLCA2\mainscatter_sigmapp10\FLAT';
fname_simul = 'LD2-26-Nov-2024_23-56-35_Final';

% variables of interest in the data
varnames = {'parsdata', 'ensdata', 'r_n_agg', 'fl'};

% load simulated data
for i = 1 : numel(varnames)
    load(strcat(fdir_simul, '\', fname_simul, '.mat'), varnames{i})
end

polyorder = 2;

%% plot dpp vs. da as a function of time

n_dat = length(r_n_agg); % number of data times to be plotted

% initialize temporal dpp vs da figure
f1 = figure(1);
f1.Position = [50, 50, 500, 600];
set(f1, 'color', 'white')

% placholders for plots & legends
plt1 = cell(n_dat+1, 1);
legtxt1 = cell(n_dat+1, 1);
legtxt1{end} = 'Olfert $\&$ Rogak (2019)';

% initialize marker colors, shapes and sizes
mc1 = colormap(hot);
cind1 = round(1 + (length(mc1) - 1) .* (0.05 : 0.7 / (n_dat - 1) : 0.75)');
mc1 = mc1(cind1,:);
ms1 = [10, 20, 15, 20, 10];
mt1 = {'^', 's', 'p', '*', 'o'};
mc1(end,:) = [236,230,61] / 255;

% generate and plot universal correlation for dpp vs. da
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
plt1{end} = plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

for i = 1 : n_dat
    
    % plot temporal data for dpp vs. da
    plt1{i} = scatter(1e9 * parsdata(i).da, 1e9 * parsdata(i).dpp,...
        ms1(i), mc1(i,:), mt1{i}, 'LineWidth', 1); 
    
    % make legends and adjust their format
    if i == 1
        legtxt1(i) = strcat('$n_\mathrm{agg}/n_\mathrm{agg_0}$ =',...
            {' '}, num2str(r_n_agg(i), '%.0f'));
    elseif ismember(i, [2,3])
        legtxt1(i) = strcat('$n_\mathrm{agg}/n_\mathrm{agg_0}$ =',...
            {' '}, num2str(r_n_agg(i), '%.1f'));
    else
        legtxt1(i) = strcat('$n_\mathrm{agg}/n_\mathrm{agg_0}$ =',...
            {' '}, num2str(r_n_agg(i), '%.2f'));
    end

    parsdata(i).r_n_agg = repelem(r_n_agg(i), length(parsdata(i).da),1);

end

dx1 = [1e9 * 0.9 * min(cat(1, parsdata.da)),...
    1e9 * 1.1 * max(cat(1, parsdata.da))];
dy1 = [1e9 * 0.9 * min(cat(1, parsdata.dpp)),...
    1e9 * 1.1 * max(cat(1, parsdata.dpp))];

% plot appearance configs
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim(dx1)
ylim(dy1)
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
box on

% generate legends
legend(cat(1, plt1{:}), legtxt1, 'interpreter', 'latex',...
    'FontSize', 11, 'orientation', 'horizontal', 'NumColumns', 2,...
    'Location', 'southoutside')

%% plot effective density as a function of time

% initialize temporal rho vs dm figure
f2 = figure(2);
f2.Position = [100, 100, 500, 600];
set(f2, 'color', 'white')

% placholders for plots & legends
plt2 = cell(n_dat+1, 1);

% initialize marker colors, shapes and sizes
mc2 = colormap(hot);
cind2 = round(1 + (length(mc2) - 1) .* (0.05 : 0.7 / (n_dat - 1) : 0.75)');
mc2 = mc2(cind2,:);
ms2 = [10, 20, 15, 20, 10];
mt2 = {'^', 's', 'p', '*', 'o'};
mc2(end,:) = [236,230,61] / 255;

% generate and plot universal correlation for rho_eff vs. dm
D_m = 2.48; % exponent
rho_eff_100 = 510; % pefactor
dm_lim_uc = [1e0 2e4];  % limits on the mobility diameter
n_dm_uc = 1e4; % number of data
uc2 = @(y) rho_eff_100 * (y / 100) .^ (D_m - 3); % on-demand function...
    % ...for the forward correlation in the mass-mobility domain...
    % ...(rho_eff in [kg/m3] as a function of dm in [nm])
r_uc2 = (dm_lim_uc(2) / dm_lim_uc(1)) ^ (1 / (n_dm_uc - 1));
dm_uc = dm_lim_uc(1) * ones(n_dm_uc,1) .* r_uc2 .^ (((1 : n_dm_uc) - 1)');
rho_eff_uc = uc2(dm_uc);
plt2{end} = plot(dm_uc, rho_eff_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

n_agg = zeros(n_dat,1); % number of aggregates in each timestep saved

rho0 = 1860; % material density

for j = 1 : n_dat
    
    n_agg = length(parsdata(j).npp);

    % initialize effective density array
    parsdata(j).rho_eff = zeros(n_agg,1);
    parsdata(j).dm = zeros(n_agg,1);
    
    % calculate mobility diameter and effective density
    for jj = 1 : n_agg
        parsdata(j).dm(jj) = TRANSP.DIAMOBIL(parsdata(j).dg(jj),...
            parsdata(j).da(jj), fl);
        parsdata(j).rho_eff(jj) = rho0 * sum(parsdata(j).pp{jj}(:,2).^3) ./...
            parsdata(j).dm(jj).^3;
    end

    % plot temporal data for rho_eff vs. dm
    plt2{j} = scatter(1e9 * parsdata(j).dm, parsdata(j).rho_eff, ms2(j),...
        mc2(j,:), mt2{j}, 'LineWidth', 1); 
    
end

dx2 = [1e9 * 0.9 * min(cat(1, parsdata.dm)),...
    1e9 * 1.1 * max(cat(1, parsdata.dm))];
dy2 = [0.9 * min(cat(1, parsdata.rho_eff)),...
    1.1 * max(cat(1, parsdata.rho_eff))];

% plot appearance configs
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim(dx2)
ylim(dy2)
xlabel('$d_\mathrm{m} [nm]$', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$\rho_\mathrm{eff} [kg/m^3]$', 'interpreter', 'latex', 'FontSize', 14)
box on

% generate legends
legend(cat(1, plt2{:}), legtxt1, 'interpreter', 'latex',...
    'FontSize', 11, 'orientation', 'horizontal', 'NumColumns', 2,...
    'Location', 'southoutside')

if ~isfolder('outputs\')
    mkdir('outputs\'); % if it doesn't exist, create the directory
end

%% plot dpp vs. da as a function of hybridity

% initialize figure 
f3 = figure(3);
f3.Position = [150, 150, 500, 600];
set(f3, 'color', 'white')

% define plotting variables
p3 = cell(7,1);
legtxt3 = cell(7,1);
legtxt3{7} = 'Olfert $\&$ Rogak (2019)';
kk3 = [1, 2, 3, 6, 11, inf];
mc3 = colormap(turbo);
ii3 = round(1 + (length(mc3) - 1) .* (0.05 : 0.9 / 4 : 0.95)');
mc3 = mc3(ii3,:);
mc3 = flip (mc3,1);
ms3 = [10, 20, 15, 20, 10];
mt3 = {'^', 's', 'p', '*', 'o'};

% plot universal correlation
p3{7,1} = plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

% concatinate the data
pars_ens.dpp = cat(1,parsdata.dpp);
pars_ens.dpp = pars_ens.dpp(:,1);
nagg_ens = length(pars_ens.dpp);
pars_ens.da = cat(1,parsdata.da);
pars_ens.nhyb = cat(1,parsdata.n_hyb);
pars_ens.pp = cat(1,parsdata.pp);

% find and remove duplicates
ij = nchoosek(1 : nagg_ens, 2);
ind_flt = zeros(nagg_ens,1);
for k = 1 : length(ij)
    if isequal(sort(unique(pars_ens.pp{ij(k,1)}(:,1))),...
            sort(unique(pars_ens.pp{ij(k,2)}(:,1))))
        ind_flt(ij(k,2)) = 1;
    end
end
ind_flt = logical(ind_flt);
pars_flt.dpp = pars_ens.dpp(~ind_flt);
pars_flt.da = pars_ens.da(~ind_flt);
pars_flt.nhyb = pars_ens.nhyb(~ind_flt);
pars_flt.pp = pars_ens.pp(~ind_flt);

% initial ensemble average pp mean diameter
pp1_ens = cell2mat(parsdata(1).pp);
dpp1_ens = geomean(pp1_ens(:,2));
sigmapp1_ens = UTILS.GEOSTD(pp1_ens(:,2));
p3{6,1} = plot(da_uc, 1e9 * repelem(dpp1_ens, length(da_uc)),...
    'Color', [0, 0, 0], 'LineStyle', ':', 'LineWidth', 2);
legtxt3(6) = strcat('$\overline{d}_\mathrm{pp,ens}$ =', {' '},...
    num2str(1e9 * dpp1_ens,'%.1f'), ' nm');

iii = cell(5,1); % index sorting placeholder based on number of internal clusters

for i = 1 : 5
    iii{i} = (pars_flt.nhyb >= kk3(i)) & (pars_flt.nhyb < kk3(i+1));
    
    % plot dpp vs. da as a function of internal cluster counts
    p3{i,1} = scatter(1e9 * pars_flt.da(iii{i}), 1e9 * pars_flt.dpp(iii{i}),...
        ms3(i), mc3(i,:), mt3{i}, 'LineWidth', 1);
    
    % make legends
    switch i
        case {1,2}
            legtxt3(i) = strcat('$n_{hyb}$ =',...
                {' '}, num2str(kk3(i), '%d'));
        case {3,4}
            legtxt3(i) = strcat(num2str(kk3(i), '%d'), {' '},...
                '$\leq n_{hyb} \leq$', {' '}, num2str(kk3(i+1)-1, '%d'));
        otherwise
            legtxt3(i) = strcat('$n_{hyb} >$', {' '},...
                num2str(kk3(i)-1, '%d'));
    end
end

box on
xlim(dx1)
ylim(dy1)
set(gca, 'FontSize', 11, 'TickLength', [0.02 0.02], 'XScale', 'log',...
    'YScale', 'log', 'TickLabelInterpreter','latex')

xlabel('$\overline{d}_\mathrm{a}$ [nm]', 'FontSize', 14, 'interpreter','latex')
ylabel('$\overline{d}_\mathrm{pp}$ [nm]', 'FontSize', 14, 'interpreter', 'latex')

legend(cat(1, p3{:})', cat(2,legtxt3(:)), 'Location',...
    'southoutside', 'FontSize', 11, 'interpreter', 'latex',...
    'NumColumns', 2)

%% plot effective density as a function of hybridity

% initialize figure 
f4 = figure(4);
f4.Position = [200, 200, 500, 600];
set(f4, 'color', 'white')

% define plotting variables
p4 = cell(6,1);
legtxt4 = cell(6,1);
legtxt4{6} = 'Olfert $\&$ Rogak (2019)';

% plot universal correlation
p4{6,1} = plot(dm_uc, rho_eff_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

% concatinate the data
pars_ens.dm = cat(1,parsdata.dm);
pars_ens.rho_eff = cat(1,parsdata.rho_eff);

% remove duplicates
pars_flt.dm = pars_ens.dm(~ind_flt);
pars_flt.rho_eff = pars_ens.rho_eff(~ind_flt);

% plot rho_eff vs. dm as a function of internal cluster counts
for i = 1 : 5
    p4{i,1} = scatter(1e9 * pars_flt.dm(iii{i}), pars_flt.rho_eff(iii{i}),...
        ms3(i), mc3(i,:), mt3{i}, 'LineWidth', 1);

    legtxt4{i} = legtxt3{i};
end

% set appearances
box on
xlim(dx2)
ylim(dy2)
set(gca, 'FontSize', 11, 'TickLength', [0.02 0.02], 'XScale', 'log',...
    'YScale', 'log', 'TickLabelInterpreter','latex')

xlabel('$d_\mathrm{m} [nm]$', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$\rho_\mathrm{eff} [kg/m^3]$', 'interpreter', 'latex', 'FontSize', 14)
legend(cat(1, p4{:})', cat(2,legtxt4(:)), 'Location',...
    'southoutside', 'FontSize', 11, 'interpreter', 'latex',...
    'NumColumns', 2)


%% fit a curve to effective density data via polynomial bayesian regression

% initialize figure 
f5 = figure(5);
f5.Position = [250, 250, 500, 500];
set(f5, 'color', 'white')

% calculate the bayesian fit based on the input order of polynomial
switch polyorder
    case 2
        [bayesfit.yfit, bayesfit.xfit, bayesfit.bounds_yfit,...
            bayesfit.afit, bayesfit.bounds_afit] =...
            UTILS.BAYESFIT_POLY2(1e9 * pars_flt.dm, pars_flt.rho_eff,...
            ones(size(pars_flt.pp)), 500);
    case 3
        [bayesfit.yfit, bayesfit.xfit, bayesfit.bounds_yfit,...
            bayesfit.afit, bayesfit.bounds_afit] =...
            UTILS.BAYESFIT_POLY3(1e9 * pars_flt.dm, pars_flt.rho_eff,...
            ones(size(pars_flt.pp)), 500);
    case 4
        [bayesfit.yfit, bayesfit.xfit, bayesfit.bounds_yfit,...
            bayesfit.afit, bayesfit.bounds_afit] =...
            UTILS.BAYESFIT_POLY4(1e9 * pars_flt.dm, pars_flt.rho_eff,...
            ones(size(pars_flt.pp)), 500);
    otherwise
        error('invalid polynomial order!!')
end

p5 = cell(3,1); % define plotting variables

% assign legends
legtxt5 = cell(3,1);
legtxt5{1} = 'Simulated aggregates';
legtxt5{2} = 'Bayesian regression';
legtxt5{3} = 'Olfert $\&$ Rogak (2019)';

% plot universal correlation
p5{3} = plot(dm_uc, rho_eff_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

p5{1} = scatter(1e9 * pars_flt.dm, pars_flt.rho_eff,...
    8, [0.35 0.35 0.35], 'o', 'LineWidth', 1);

% mean slope
p5{2} = plot(bayesfit.xfit, bayesfit.yfit, 'Color', hex2rgb('#659287'),...
    'LineWidth', 2);
% shaded 95% CI
fill([bayesfit.xfit; flipud(bayesfit.xfit)], [bayesfit.bounds_yfit(:,1);...
      flipud(bayesfit.bounds_yfit(:,2))], hex2rgb('#659287'),...
      'EdgeColor', 'none', 'FaceAlpha', 0.3);

% set appreances
box on
xlim(dx2)
ylim(dy2)
set(gca, 'FontSize', 11, 'TickLength', [0.02 0.02], 'XScale', 'log',...
    'YScale', 'log', 'TickLabelInterpreter','latex')
xlabel('$d_\mathrm{m} [nm]$', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$\rho_\mathrm{eff} [kg/m^3]$', 'interpreter', 'latex', 'FontSize', 14)
legend(cat(1, p5{:})', cat(2,legtxt5(:)), 'Location',...
    'southoutside', 'FontSize', 11, 'interpreter', 'latex',...
    'NumColumns', 2)

%% plot mass-mobility exponent (slope of effective density in log-log scale)

% initialize figure 
f6 = figure(6);
f6.Position = [300, 300, 500, 500];
set(f6, 'color', 'white')

p6 = cell(3,1); % define plotting variables

% assign legends
legtxt6 = cell(3,1);
legtxt6{1} = 'Regressed simulation';
legtxt6{2} = 'Olfert $\&$ Rogak (2019)';
legtxt6{3} = 'DLCA limit';

% mean slope
p6{1} = plot(bayesfit.xfit, bayesfit.afit, 'Color', hex2rgb('#659287'),...
    'LineWidth', 2);
hold on
% shaded 95% CI
fill([bayesfit.xfit; flipud(bayesfit.xfit)], [bayesfit.bounds_afit(:,1);...
      flipud(bayesfit.bounds_afit(:,2))], hex2rgb('#659287'),...
      'EdgeColor', 'none', 'FaceAlpha', 0.3);

% plot exponent of universal correlation
p6{2} = plot(bayesfit.xfit, 2.48 * ones(size(bayesfit.xfit)),...
    'Color', [0.4940 0.1840 0.5560], 'LineStyle', '-.', 'LineWidth', 3);

% plot exponent limit for diffusion-limited cluster aggregation
p6{3} = plot(bayesfit.xfit, 1.78 * ones(size(bayesfit.xfit)),...
    'Color', [0, 0, 0], 'LineStyle', ':', 'LineWidth', 3);

xlim([min(bayesfit.xfit), 1000])
ylim([1.5, 2.8])
set(gca, 'FontSize', 11, 'TickLength', [0.02 0.02], 'XScale', 'log',...
    'YScale', 'log', 'TickLabelInterpreter','latex')
xlabel('$d_\mathrm{m} [nm]$', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$D_\mathrm{m} [-]$', 'interpreter', 'latex',...
    'FontSize', 14)
legend(cat(1, p6{:})', cat(2,legtxt6(:)), 'Location',...
    'southoutside', 'FontSize', 11, 'interpreter', 'latex',...
    'NumColumns', 2)

%% save plots

% make a directory to save outputs
dir_out = strcat('outputs\', 'postLD2-', datestr(datetime('now')),...
    '_', fname_simul, '\');
dir_out = regexprep(dir_out, ':', '-');
dir_out = regexprep(dir_out, ' ', '_');
if ~isfolder(dir_out)
    mkdir(dir_out); % if it doesn't exist, create the directory
end

% save worksapce
save(strcat(dir_out, 'Post_', fname_simul, '.mat'))

% print figures
exportgraphics(f1, strcat(dir_out, 'dpp-da-temporal.jpg'),...
    'BackgroundColor','none', 'Resolution', 300)
exportgraphics(f2, strcat(dir_out, 'rho-dm-temporal.jpg'),...
    'BackgroundColor','none', 'Resolution', 300)
exportgraphics(f3, strcat(dir_out, 'dpp-da-hybrid.jpg'),...
    'BackgroundColor','none', 'Resolution', 300)
exportgraphics(f4, strcat(dir_out, 'rho-dm-hybrid.jpg'),...
    'BackgroundColor','none', 'Resolution', 300)
exportgraphics(f6, strcat(dir_out, 'Dm-dm.jpg'),...
    'BackgroundColor','none', 'Resolution', 300)
