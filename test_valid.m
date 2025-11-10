clear; clc; close all;

% load scaled uniform aggregates
fdir_dpp_da = 'D:\Hamed\CND\PhD\Weekly\2025\Week_of_27OCT25';
fname_dpp_da = 'BiggerRepo-Scat135-Real1-Frac0055-27OCT25';
vars_dpp_da = {'pars_out'};
load(fullfile(fdir_dpp_da, fname_dpp_da), vars_dpp_da{:});

% coefficients for Olfert and Rogak (2019)'s universal correlation
D_TEM = 0.35; % exponent
dpp100 = 17.8; % pefactor

% load aggregate sizes and convert units
da = 1e9 * pars_out.da; % projected area diameter
dpp = 1e9 * pars_out.dpp(:,1); % primary particle diameter

% get 3 evenly log-spaced targets for da
log_da = log10(da);
log_da_targets = linspace(min(log_da), max(log_da), 5);
log_da_targets = log_da_targets(2:4);

%%% find indices of closest and furthest points to universal correlation

% get 3 log-spaced da values
log_da = log10(da);
log_da_targets = linspace(min(log_da), max(log_da), 5);
log_da_targets = log_da_targets(2:4);


id = zeros(3, 3); % indices of closest and furthest points to targets

for i = 1:3 % for each da target

    % find points within da window
    da_window = 0.05;  % adjust as needed (unit: nm)
    nearby = find(abs(log_da - log_da_targets(i)) < da_window);
    
    % distance from line for these points
    dist = log10(dpp(nearby)) - (D_TEM*log10(da(nearby)/100) + log10(dpp100));
    
    % closest to line
    [~, close_id] = min(abs(dist));
    id(i,1) = nearby(close_id);
    
    % furthest above line
    [~, top_id] = max(dist);
    id(i,2) = nearby(top_id);
    
    % furthest below line
    [~, bot_id] = min(dist);
    id(i,3) = nearby(bot_id);

end

%%% draw selected points in dpp vs da space against universal correlation

% initialize figure
f1 = figure(1);
f1.Position = [100, 100, 450, 450];
set(f1, 'color', 'white');

% draw universal correlation
n_uc = 100;
da_lim_uc = [min(da) max(da)];
r_uc = (da_lim_uc(2) / da_lim_uc(1)) ^ (1 / (n_uc - 1));
da_uc = min(da) * ones(n_uc,1) .* r_uc .^ (((1 : n_uc) - 1)');
dpp_uc = dpp100 * (da_uc/100) .^ D_TEM;
plt_uc = plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

% appearance configs for dpp vs npp subplot
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([min(da) max(da)])
ylim([min(dpp)/1.1 max(dpp)*1.1])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 16)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 16)
box on

% draw points found in the selection loop
colors = [0.7 0.9 0.7; 0.9 0.7 0.7; 0.7 0.7 0.9]; % pastel green, red, blue
markers = ['o', 's', '^'];
sizes = [50, 60, 40];
for i = 1:3
    for j = 1:3
        scatter(da(id(i,j)), dpp(id(i,j)), sizes(i), colors(j,:), markers(i), ...
            'LineWidth', 2);
    end
end

% legend for universal correlation
legend(plt_uc, 'Olfert $\&$ Rogak (2019)', 'interpreter', 'latex',...
    'FontSize', 12, 'Location', 'northoutside');

%%% map selected points to effective density vs. mobility diameter domain

% initialize figure
f2 = figure(2);
f2.Position = [200, 200, 450, 450];
set(f2, 'color', 'white');

% address of second-stage langevin dynamics data to be imported
fdir_LD2 = 'D:\Hamed\CND\PhD\Weekly\2025\Week_of_27OCT25\Results\LD2-28-Oct-2025_BiggerRepo-Scat135-Real1-Frac0055-27OCT25';
fname_LD2 = 'LD2_30-Oct-2025_12-50-54_Final';
varnames_LD2 = {'fl'};
load(fullfile(fdir_LD2, fname_LD2), varnames_LD2{:});

% calculate mobility diameter for all aggregates
dm = TRANSP.DIAMOBIL(pars_out.dg, pars_out.da, fl);

% calculate their effective density
rho_eff = cellfun(@(x) 1860 * (pi/6) * sum(x(:,2).^3), pars_out.pp) ./ dm.^3;

dm = 1e9*dm; % convert to [nm]

% draw universal correlation (for effective density)
D_m = 2.48; % exponent
rho_eff_100 = 510; % pefactor
dm_lim_uc = [min(dm) max(dm)];
r_uc_2 = (dm_lim_uc(2) / dm_lim_uc(1)) ^ (1 / (n_uc - 1));
dm_uc = min(dm) * ones(n_uc,1) .* r_uc_2 .^ (((1 : n_uc) - 1)');
rho_eff_uc = rho_eff_100 * (dm_uc/100) .^ (D_m-3);
plt_uc_2 = plot(dm_uc, rho_eff_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

% appearance configs for dpp vs npp subplot
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([min(dm) max(dm)])
ylim([min(rho_eff)/1.1 max(rho_eff)*1.1])
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 16)
ylabel('$\rho_\mathrm{eff}$ [nm]', 'interpreter', 'latex', 'FontSize', 16)
box on

% draw points found in the selection loop
for i = 1:3
    for j = 1:3
        scatter(dm(id(i,j)), rho_eff(id(i,j)), sizes(i), colors(j,:), markers(i), ...
            'LineWidth', 2);
    end
end

% legend
legend(plt_uc_2, 'Olfert $\&$ Rogak (2019)', 'interpreter', 'latex',...
    'FontSize', 12, 'Location', 'northoutside');

