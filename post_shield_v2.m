%% Post-process shielding bias for TWO datasets (side-by-side comparison)
% Figures:
%   Fig.1  : S_pp distributions (KDE-only), two panels A|B (left y-label only)
%   Fig.2  : <S_pp> vs n_pp (scatter), two panels A|B with ONE shared legend
%   Fig.3  : Detailed renders merged into ONE A4-friendly figure (6x3):
%            rows 1–3 = case A (observable | all | screened),
%            rows 4–6 = case B; titles = case labels only; original titles
%            restored as SUBTITLES; in rows 1 and 4, the *middle* panel
%            title is just the case label ("lacey"/"compact").
%   Fig.4  : 2x2 (A|B columns). Top: boxplots (spp* = 0.5) with ensemble
%            baseline line. Bottom: bias curves via KDE (shared support +
%            bandwidth). Left y-labels only.
%   Fig.5  : 3x2. Row1: dpp(2D)/dpp(3D) vs d_a (scatter). Row2:
%            dpp(2D)/dpp(3D) vs n_agg/(n_agg)_2 (box). Row3:
%            sigma_pp(2D)/sigma_pp(3D) vs n_agg/(n_agg)_2 (box).
%            Horizontal y=1 guide lines and gentle bands (for clarity).
%
% Global rules you requested:
% - Titles in ALL figures are only the case label (“lacey”/“compact”),
%   with original per-panel text moved to SUBTITLES where you asked.
% - For A|B figures, print the y-label only on the LEFT panel.
% - Clean, non-overlapping mean labels in Fig.1; KDE-only (no hist).
%
% Performance:
% - Vectorized concatenations (once per snapshot).
% - KDEs use shared support and bandwidth per panel for stability.
% - No expensive styling in tight loops (kept where worth it visually).

clc
clear
clf('reset')
close all
warning('off')

%% ----------------------- User-configurable inputs -----------------------

% --- Case labels (for titles/legends) ---
case_label_A = 'lacey';    % default per your request
case_label_B = 'compact';

% --- Input files (CASE A) ---
fdir_in_A  = 'D:\Hamed\CND\PhD\Publication\DLCA2\outputs\Shield_2025-08-09_18-03-51';
fname_in_A = 'Shield_2025-08-09_18-03-51';

% --- Input files (CASE B) ---
% (Set these to your second dataset)
fdir_in_B  = 'D:\Hamed\CND\PhD\Publication\DLCA2\outputs\Shield_2025-08-09_18-03-51';
fname_in_B = 'Shield_2025-08-09_18-03-51';

% --- Variables to import ---
varnames = {'parsdata'}; % must contain `parsdata` as in your schema

% --- Shielding thresholds to evaluate ---
spp_star = [1, 0.75, 0.5];
i_spp_star = 3;  % index selecting spp_star for some plots (e.g., bias)

% --- Appearance options ---
font_base   = 16;
font_big    = 18;
font_title  = 16;  % panel title size
font_subttl = 12;  % subtitle size (was original per-panel title)
ms_vec      = [8, 16, 12, 16, 8];   % marker sizes cycling
mt_vec      = {'^','s','p','*','o'}; % marker types cycling

% --- Fig.1 (spp distributions) options ---
fig1_use_kde = true;        % KDE only (per your request, hist removed)
fig1_mean_font = 14;        % mean label font (reduced to avoid overlap)

% --- Fig.4 bias (bottom row) options ---
bias_use_kde = true;        % use KDE for bias curves (less noisy)
bias_n_grid  = 256;         % log-grid points for KDE support

% --- Gentle band for ratio plots (Fig.5 rows 2 & 3) ---
ratio_band1 = 0.02;  % ±2% band around y=1 (light shading)

%% --------------------------- Load both cases ----------------------------

% Load CASE A
load(fullfile(fdir_in_A, strcat(fname_in_A, '.mat')), varnames{1});
parsA = parsdata; clear parsdata

% Load CASE B
load(fullfile(fdir_in_B, strcat(fname_in_B, '.mat')), varnames{1});
parsB = parsdata; clear parsdata

%% ----------------------- Common snapshot bookkeeping --------------------

n_shot_A = numel(parsA);
n_shot_B = numel(parsB);

% Generate x tick labels per case (independent if counts differ)
xlbl_A = gen_xticklabels(parsA);
xlbl_B = gen_xticklabels(parsB);

%% ----------------------- Colors & markers (shared) ----------------------

clr_all = colormap(hot);
cind_A = round(1 + (length(clr_all) - 1) .* (0.05 : 0.7 / (max(n_shot_A,5) - 1) : 0.75)');
clrA   = clr_all(cind_A(1:n_shot_A),:); if ~isempty(clrA), clrA(end,:) = [236,230,61]/255; end

cind_B = round(1 + (length(clr_all) - 1) .* (0.05 : 0.7 / (max(n_shot_B,5) - 1) : 0.75)');
clrB   = clr_all(cind_B(1:n_shot_B),:); if ~isempty(clrB), clrB(end,:) = [236,230,61]/255; end

msA = ms_vec(1:min(n_shot_A,numel(ms_vec)));
mtA = mt_vec(1:min(n_shot_A,numel(mt_vec)));
msB = ms_vec(1:min(n_shot_B,numel(ms_vec)));
mtB = mt_vec(1:min(n_shot_B,numel(mt_vec)));

%% ---------------- Vectorized precomputation per case/snapshot -----------

SA = prep_case(parsA);
SB = prep_case(parsB);

%% ----------------------- Fig.1: spp distributions -----------------------
% Two panels: left = Case A, right = Case B
% KDE-only, with clean mean labels and no overlap. Left y-label only.

f1 = figure(1); set(f1, 'Color','w','Position',[100 100 1100 520]);
tl1 = tiledlayout(1,2,"TileSpacing","compact","Padding","compact");

% Left tile (Case A)
axA = nexttile(tl1,1);
plot_spp_panel_kde(axA, SA, clrA, xlbl_A, case_label_A, font_base, fig1_mean_font);

% Right tile (Case B)
axB = nexttile(tl1,2);
plot_spp_panel_kde(axB, SB, clrB, xlbl_B, case_label_B, font_base, fig1_mean_font);

% Left y-label only
ylabel(axA, '$S_\mathrm{pp}^\mathrm{(i)}$ [-]', 'Interpreter','latex', 'FontSize', font_big);
ylabel(axB, ''); % suppress duplicate
% x-limits: [-0.4, n_shot+0.4], tick labels already set in the panels

%% --------------- Fig.2: <S_pp> vs npp (scatter), A | B ------------------

f2 = figure(2); set(f2,'Color','w','Position',[140 140 1100 520]);
tl2 = tiledlayout(1,2,"TileSpacing","compact","Padding","compact");

[scatsA, legtxtA] = plot_sagg_vs_npp(nexttile(tl2,1), SA, parsA, clrA, msA, mtA, xlbl_A, case_label_A, font_base);
[scatsB, legtxtB] = plot_sagg_vs_npp(nexttile(tl2,2), SB, parsB, clrB, msB, mtB, xlbl_B, case_label_B, font_base);

% Left y-label only
ylabel(nexttile(tl2,1), '$S_\mathrm{pp}$ [-]', 'Interpreter','latex','FontSize', font_big);
ylabel(nexttile(tl2,2), '');

% Shared legend (one only). Prefer outside south of layout.
% Use the left panel's legend strings (they mirror the snapshot labels).
% keep the handles of the two axes when you plot
ax2L = nexttile(tl2,1);
[scatsA, legtxtA] = plot_sagg_vs_npp(ax2L, SA, parsA, clrA, msA, mtA, xlbl_A, case_label_A, font_base);

ax2R = nexttile(tl2,2);
[scatsB, legtxtB] = plot_sagg_vs_npp(ax2R, SB, parsB, clrB, msB, mtB, xlbl_B, case_label_B, font_base);

% one shared legend, anchored to an axes, then positioned in the layout
lgd = legend(ax2L, scatsA, legtxtA, 'Interpreter','latex','FontSize',font_base, ...
    'NumColumns', min(3, numel(scatsA)));
lgd.Layout.Tile = 'south';

% y-label only on the left
ylabel(ax2L, '$S_\mathrm{pp}$ [-]', 'Interpreter','latex','FontSize', font_big);
ylabel(ax2R, '');

%% ---------------- Fig.3: detailed renders merged (6x3) ------------------
% Rows 1–3: case A; Rows 4–6: case B
% Columns: Observable | All | Screened
% Titles = case labels only; restore your original per-panel "titles" as subtitles.
% In rows 1 and 4, the middle panel title is just the case label.

% Selection of snapshots/aggregates to render
ii1 = [1, 2, 3]; % snapshots to render (as before)
% Find example aggregates per your matcher (guard for each case)
jjA = UTILS.MATCHAGG(vertcat(parsA(ii1)), 'npp','sigmapp','n_hyb', 0.8,80,50.0,10000,38,300);
jjB = UTILS.MATCHAGG(vertcat(parsB(ii1)), 'npp','sigmapp','n_hyb', 0.8,80,50.0,10000,38,300);

f3 = figure(3); set(f3,'Color','w','Position',[180 60 1000 1350]); % A4-friendly portrait
tl3 = tiledlayout(6,3,"TileSpacing","compact","Padding","compact");

% CASE A rows 1–3
for t = 1:numel(ii1)
    i = ii1(t);
    % Observable (spp < 0.5)
    nexttile(tl3, 3*(t-1) + 1);
    render_pp_subset(parsA(i).pp{jjA(t)}, parsA(i).spp{jjA(t)}, [0,0.5], 'blue');
    if t==1
        title(case_label_A, 'Interpreter','latex', 'FontSize', font_title);
        subtitle('Observable primaries (spp < 0.5)', 'Interpreter','latex', 'FontSize', font_subttl);
    else
        title(case_label_A, 'Interpreter','latex', 'FontSize', font_title);
        subtitle('Observable primaries (spp < 0.5)', 'Interpreter','latex', 'FontSize', font_subttl);
    end

    % All
    nexttile(tl3, 3*(t-1) + 2);
    render_pp_subset(parsA(i).pp{jjA(t)}, parsA(i).spp{jjA(t)}, [0,1], 'purple');
    if t==1
        title(case_label_A, 'Interpreter','latex', 'FontSize', font_title); % (row 1, middle) title = case label
        subtitle('All primaries', 'Interpreter','latex', 'FontSize', font_subttl);
    else
        title(case_label_A, 'Interpreter','latex', 'FontSize', font_title);
        subtitle('All primaries', 'Interpreter','latex', 'FontSize', font_subttl);
    end

    % Screened (spp >= 0.5)
    nexttile(tl3, 3*(t-1) + 3);
    render_pp_subset(parsA(i).pp{jjA(t)}, parsA(i).spp{jjA(t)}, [0.5,1], 'orange');
    title(case_label_A, 'Interpreter','latex', 'FontSize', font_title);
    subtitle('Screened primaries (spp \ge 0.5)', 'Interpreter','latex', 'FontSize', font_subttl);
end

% CASE B rows 4–6
offset = 9; % 3 rows * 3 cols = 9 tiles above
for t = 1:numel(ii1)
    i = ii1(t);
    % Observable
    nexttile(tl3, offset + 3*(t-1) + 1);
    render_pp_subset(parsB(i).pp{jjB(t)}, parsB(i).spp{jjB(t)}, [0,0.5], 'blue');
    title(case_label_B, 'Interpreter','latex', 'FontSize', font_title);
    subtitle('Observable primaries (spp < 0.5)', 'Interpreter','latex', 'FontSize', font_subttl);

    % All (row 4 middle panel: title = case label)
    nexttile(tl3, offset + 3*(t-1) + 2);
    render_pp_subset(parsB(i).pp{jjB(t)}, parsB(i).spp{jjB(t)}, [0,1], 'purple');
    if t==1
        title(case_label_B, 'Interpreter','latex', 'FontSize', font_title); % row 4 middle
    else
        title(case_label_B, 'Interpreter','latex', 'FontSize', font_title);
    end
    subtitle('All primaries', 'Interpreter','latex', 'FontSize', font_subttl);

    % Screened
    nexttile(tl3, offset + 3*(t-1) + 3);
    render_pp_subset(parsB(i).pp{jjB(t)}, parsB(i).spp{jjB(t)}, [0.5,1], 'orange');
    title(case_label_B, 'Interpreter','latex', 'FontSize', font_title);
    subtitle('Screened primaries (spp \ge 0.5)', 'Interpreter','latex', 'FontSize', font_subttl);
end

%% ---------- Compute dpp_2d (by thresholds) + per-agg stats for A/B ------

n_spp_star = numel(spp_star);

[dpp_2d_A, sig_2d_A, n_dpp2d_A] = compute_obs_stats(parsA, SA, spp_star);
[dpp_2d_B, sig_2d_B, n_dpp2d_B] = compute_obs_stats(parsB, SB, spp_star);

% Per-agg 3D sigma (for Fig.5 third row)
sig3D_agg_A = peragg_sig3D(parsA);
sig3D_agg_B = peragg_sig3D(parsB);

%% ---------------- Fig.4: Bias figure duplicated A/B (2×2) ---------------
% Top: boxplots (spp* = 0.5) + ensemble baseline line
% Bottom: bias curves (KDE with shared support/bandwidth). Left y-labels only.

f4 = figure(4); set(f4,'Color','w','Position',[260 0 1200 1250]);
tl4 = tiledlayout(2,2,"TileSpacing","compact","Padding","compact");

% A: boxplots
ax41 = nexttile(tl4,1);
plot_box_tile(ax41, dpp_2d_A, SA.n_shot, xlbl_A, SA.dpp3D_ens, clrA, font_base, case_label_A);
ylabel(ax41, '$d_\mathrm{pp}^{(i)}$ [nm]', 'Interpreter','latex','FontSize', font_big);

% B: boxplots
ax42 = nexttile(tl4,2);
plot_box_tile(ax42, dpp_2d_B, SB.n_shot, xlbl_B, SB.dpp3D_ens, clrB, font_base, case_label_B);
ylabel(ax42, ''); % right tile, no y-label

% A: bias curves (KDE, shared support + bandwidth)
ax43 = nexttile(tl4,3);
plot_bias_kde(ax43, dpp_2d_A, SA.n_shot, clrA, font_base, bias_n_grid);
ylabel(ax43, '$f^{(2D)}/f^{(3D)}$', 'Interpreter','latex','FontSize', font_big);

% B: bias curves
ax44 = nexttile(tl4,4);
plot_bias_kde(ax44, dpp_2d_B, SB.n_shot, clrB, font_base, bias_n_grid);
ylabel(ax44, '');

%% ---------------- Fig.5: 3×2 comparison (from old Fig.7) ----------------
% Row 1: dpp(2D)/dpp(3D) vs d_a   (A | B)
% Row 2: dpp(2D)/dpp(3D) vs n_agg/(n_agg)_2   (A | B)
% Row 3: sigma_pp(2D)/sigma_pp(3D) vs n_agg/(n_agg)_2  (A | B)
% Add y=1 horizontal guide and gentle ±band for clarity.

f5 = figure(5); set(f5,'Color','w','Position',[320 40 1200 1000]);
tl5 = tiledlayout(3,2,"TileSpacing","compact","Padding","compact");

% Row 1
ax51 = nexttile(tl5,1);
plot_ratio_vs_da(ax51, parsA, dpp_2d_A, clrA, case_label_A, font_base);
ylabel(ax51, '$d_\mathrm{pp}^{(2D)} / d_\mathrm{pp}^{(3D)}$ [-]', 'Interpreter','latex','FontSize', font_big);

ax52 = nexttile(tl5,2);
plot_ratio_vs_da(ax52, parsB, dpp_2d_B, clrB, case_label_B, font_base);
ylabel(ax52, '');

% Row 2
ax53 = nexttile(tl5,3);
plot_ratio_vs_snap(ax53, parsA, dpp_2d_A, clrA, xlbl_A, case_label_A, font_base, ratio_band1);
ylabel(ax53, '$d_\mathrm{pp}^{(2D)} / d_\mathrm{pp}^{(3D)}$ [-]', 'Interpreter','latex','FontSize', font_big);

ax54 = nexttile(tl5,4);
plot_ratio_vs_snap(ax54, parsB, dpp_2d_B, clrB, xlbl_B, case_label_B, font_base, ratio_band1);
ylabel(ax54, '');

% Row 3
ax55 = nexttile(tl5,5);
plot_sig_ratio_vs_snap(ax55, parsA, sig_2d_A, sig3D_agg_A, clrA, xlbl_A, case_label_A, font_base, ratio_band1);
ylabel(ax55, '$\sigma_\mathrm{pp}^{(2D)} / \sigma_\mathrm{pp}^{(3D)}$ [-]', 'Interpreter','latex','FontSize', font_big);

ax56 = nexttile(tl5,6);
plot_sig_ratio_vs_snap(ax56, parsB, sig_2d_B, sig3D_agg_B, clrB, xlbl_B, case_label_B, font_base, ratio_band1);
ylabel(ax56, '');

%% ------------------------ Save outputs & workspace -----------------------

% Output dir
dir0_out = datestr(datetime('now'));
dir0_out = regexprep(dir0_out, ':', '-');
dir0_out = regexprep(dir0_out, ' ', '_');
dir_out = fullfile('outputs', ['PostShield2Case_', dir0_out, filesep]);
if ~isfolder(dir_out), mkdir(dir_out); end

save(fullfile(dir_out, ['PostShield2Case_', dir0_out, '.mat']), ...
     'parsA','parsB','SA','SB','dpp_2d_A','dpp_2d_B','sig_2d_A','sig_2d_B');

exportgraphics(f1, fullfile(dir_out,'Fig1_spp_distributions.jpg'), 'BackgroundColor','none','Resolution',300);
exportgraphics(f2, fullfile(dir_out,'Fig2_sagg_vs_npp.jpg'),      'BackgroundColor','none','Resolution',300);
exportgraphics(f3, fullfile(dir_out,'Fig3_renders_6x3.jpg'),      'BackgroundColor','none','Resolution',300);
exportgraphics(f4, fullfile(dir_out,'Fig4_bias_comp_2x2.jpg'),    'BackgroundColor','none','Resolution',300);
exportgraphics(f5, fullfile(dir_out,'Fig5_ratios_3x2.jpg'),       'BackgroundColor','none','Resolution',300);

%% ============================= Local functions ==========================

function xlbl = gen_xticklabels(pars)
    n_shot = numel(pars);
    xlbl = cell(n_shot,1);
    for i = 1:n_shot
        r = pars(i).r_n_agg(1);
        if i == 1
            xlbl{i} = sprintf('%.0f', r);
        elseif ismember(i,[2,3])
            xlbl{i} = sprintf('%.1f', r);
        else
            xlbl{i} = sprintf('%.2f', r);
        end
    end
end

function S = prep_case(pars)
    S.n_shot = numel(pars);
    S.pp_all = cell(S.n_shot,1);
    S.spp_all = cell(S.n_shot,1);
    S.nagg = zeros(S.n_shot,1);
    S.sagg = cell(S.n_shot,1);
    for i = 1:S.n_shot
        S.nagg(i)    = numel(pars(i).pp);
        S.pp_all{i}  = cat(1, pars(i).pp{:});   % columns: uses col 2 = dpp
        S.spp_all{i} = cat(1, pars(i).spp{:});
        tmp = zeros(S.nagg(i),1);
        for k = 1:S.nagg(i), tmp(k) = mean(pars(i).spp{k}); end
        S.sagg{i} = tmp;
    end
    all_pp = cell2mat(S.pp_all);
    S.dpp3D_ens    = geomean(all_pp(:,2), 'omitnan');
    S.sigmapp3D_ens = geostd_local(all_pp(:,2));
end

function gsd = geostd_local(x)
    x = x(~isnan(x) & x>0);
    if isempty(x), gsd = NaN; return; end
    lx = log(x);
    gsd = exp(std(lx, 0, 'omitnan'));
end

function plot_spp_panel_kde(ax, S, clr, xlbl, label_case, font_base, mean_font)
    axes(ax); cla; hold on
    n_shot = S.n_shot;
    xi = linspace(0,1,512);

    % Precompute KDEs and means
    dens = zeros(n_shot, numel(xi));
    mu_spp = zeros(n_shot,1);
    maxf = 0;
    for i = 1:n_shot
        s = S.spp_all{i};
        mu_spp(i) = mean(s,'omitnan');
        dens(i,:) = ksdensity(s, xi, 'Function','pdf');
        maxf = max(maxf, max(dens(i,:)));
    end

    % Automatic vertical scaling to avoid overlap
    if maxf > 0
        scale = -0.85 / maxf; % pack curves well within row gap
    else
        scale = -0.25;
    end

    % Draw reference (snapshot 1) lightly under each row; then each curve
    for i = 1:n_shot
        % Light reference of snap 1 for comparison
        plot(scale*dens(1,:) + (i-0.5), xi, 'Color', [clr(1,:) 0.35], 'LineWidth', 0.5);
        fill([scale*dens(1,:) + (i-0.5), (i-0.5)*ones(size(xi))], [xi, fliplr(xi)], ...
             clr(1,:), 'FaceAlpha', 0.08, 'EdgeColor', 'none');

        % KDE for current snapshot
        plot(scale*dens(i,:) + (i-0.5), xi, 'Color', clr(i,:), 'LineWidth', 2);
        fill([scale*dens(i,:) + (i-0.5), (i-0.5)*ones(size(xi))], [xi, fliplr(xi)], ...
             clr(i,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');

        % Tick placeholder and mean labels
        scatter(i-0.5, 1.005, 1, 'w', 'filled'); % invisible anchor
        if i == 1
            text((i - 0.82), 1.035, sprintf('$\\langle S_\\mathrm{pp}\\rangle$ = %.2f', mu_spp(i)), ...
                'Interpreter','latex','HorizontalAlignment','center','FontSize', mean_font);
        else
            text((i - 0.48), 1.035, sprintf('%.2f', mu_spp(i)), ...
                'Interpreter','latex','HorizontalAlignment','center','FontSize', mean_font);
        end
    end

    % Axes cosmetics
    xticks((1:n_shot)-0.5); xticklabels(xlbl);
    xlim([-0.4, n_shot+0.4]); ylim([0, 1]);
    box on
    set(gca,'TickLabelInterpreter','latex','FontSize',font_base,'TickLength',[0.02 0.02]);
    xlabel('$n_\mathrm{agg}/(n_\mathrm{agg})_2$ [-]','Interpreter','latex','FontSize', font_base+4);
    % Title = case label only
    title(label_case, 'Interpreter','latex','FontSize', font_base+2);
end

function [scats, legtxt] = plot_sagg_vs_npp(ax, S, pars, clr, ms, mt, xlbl, label_case, font_base)
    axes(ax); cla; hold on
    n_shot = S.n_shot;
    scats = gobjects(n_shot,1);
    legtxt = cell(n_shot,1);

    for i = 1:n_shot
        legtxt{i} = sprintf('$n_\\mathrm{agg}/(n_\\mathrm{agg})_2$ = %s', xlbl{i});
        sagg_i = S.sagg{i};
        npp_i  = cellfun(@(x) size(x,1), pars(i).pp(:));
        scats(i) = scatter(npp_i, sagg_i, ms(min(i,numel(ms))), ...
            clr(i,:), mt{min(i,numel(mt))}, 'LineWidth', 1);
    end

    % bounds
    npp_all  = cat(1, pars.npp);
    sagg_all = cat(1, S.sagg{:});
    xlim([0.7*min(npp_all), 1.3*max(npp_all)]);
    ylim([0.95*min(sagg_all), 1.05*max(sagg_all)]);

    % cosmetics
    set(gca,'XScale','log','YScale','log','TickLabelInterpreter','latex', ...
        'FontSize',font_base,'TickLength',[0.02 0.02])
    xlabel('$n_\mathrm{pp}$ [-]','Interpreter','latex','FontSize',font_base+4)
    % y-label applied only on left tile in caller
    yticks([linspace(0.01,0.1,10), linspace(0.2,1,9)])
    box on
    title(label_case, 'Interpreter','latex','FontSize', font_base+2);
end

function render_pp_subset(pp, spp, clim, whichmap)
    % Inputs:
    %  pp  : N x 5 (dpp in col 2; positions in cols 3:5)
    %  spp : N x 1 shielding
    %  clim: [smin smax] color limits
    %  whichmap: 'blue' | 'purple' | 'orange'
    if strcmpi(whichmap,'blue')
        cm = UTILS.CUSTOMABYSSMAP('blue');
    elseif strcmpi(whichmap,'purple')
        cm = UTILS.CUSTOMABYSSMAP('purple');
    else
        cm = UTILS.CUSTOMABYSSMAP('orange');
    end
    cm = flip(cm,1);
    kk = (spp >= clim(1)) & (spp <= clim(2));
    UTILS.PLOTPP_CONTINUOUS(pp(kk,3), pp(kk,4), pp(kk,5), pp(kk,2), spp(kk), ...
        struct('cc','on','cm',cm,'clim',clim,'ft',0.9));
    axis vis3d
end

function [dpp_2d, sig_2d, n_dpp2d] = compute_obs_stats(pars, S, spp_star)
    n_spp = numel(spp_star);
    dpp_2d = cell(S.n_shot, n_spp, 2);
    sig_2d  = cell(S.n_shot, n_spp);
    n_dpp2d = zeros(S.n_shot, n_spp);

    for i = 1:S.n_shot
        for j = 1:n_spp
            dpp_2d{i,j,2} = zeros(S.nagg(i),1);
            sig_2d{i,j}   = zeros(S.nagg(i),1);
            mask_all = cell(S.nagg(i),1);
            for k = 1:S.nagg(i)
                mask_all{k} = pars(i).spp{k} <= spp_star(j);
                x = pars(i).pp{k}(mask_all{k},2);
                dpp_2d{i,j,2}(k) = geomean(x,'omitnan');
                sig_2d{i,j}(k)   = geostd_local(x);
            end
            dpp_all = cat(1, pars(i).pp{:});
            dpp_2d{i,j,1} = dpp_all(cat(1, mask_all{:}), 2);
            n_dpp2d(i,j)  = numel(dpp_2d{i,j,1});
        end
    end
end

function sig3D_agg = peragg_sig3D(pars)
    n_shot = numel(pars);
    sig3D_agg = cell(n_shot,1);
    for i = 1:n_shot
        na = numel(pars(i).pp);
        s = zeros(na,1);
        for k = 1:na
            s(k) = geostd_local(pars(i).pp{k}(:,2));
        end
        sig3D_agg{i} = s;
    end
end

function plot_box_tile(ax, dpp_2d_case, n_shot, xlbl, dpp3D_ens, clr, font_base, case_label)
    axes(ax); cla; hold on
    for i = 1:n_shot
        d_obs = 1e9 * dpp_2d_case{i,3,1}; % spp* = 0.5 (index 3)
        if isempty(d_obs) || all(isnan(d_obs)), continue; end
        b = boxplot(d_obs, 'Positions', i, 'Notch','on', 'Symbol','o', 'Widths', 0.5);
        boxObj = findobj(b,'Tag','Box');
        patch(get(boxObj,'XData'), get(boxObj,'YData'), clr(i,:), 'FaceAlpha',0.3, ...
              'EdgeColor', clr(i,:), 'LineWidth', 1);
        set(findobj(b,'Tag','Median'), 'Color', clr(i,:), 'LineWidth', 2);
        set(findobj(b,'Tag','Outliers'), 'MarkerEdgeColor', clr(i,:), 'MarkerSize', 3);
        set(findobj(b,'Tag','Upper Whisker'), 'LineStyle','-');
        set(findobj(b,'Tag','Lower Whisker'), 'LineStyle','-');
    end
    % Ensemble 3D mean guide (your earlier top-row guide)
    plot(linspace(0, n_shot+1, 100), 1e9*repmat(dpp3D_ens,1,100), ...
         'k:', 'LineWidth',1.5);
    xticks(1:n_shot); xticklabels(xlbl);
    set(gca,'YScale','log','TickLabelInterpreter','latex','FontSize',font_base);
    xlim([0.25, n_shot+0.75]);
    yticks([1 2 3 4 5 6 8 10 20 30 40 50 60 80 100]);
    xlabel('$n_\mathrm{agg}/(n_\mathrm{agg})_2$ [-]','Interpreter','latex','FontSize',font_base+4);
    title(case_label, 'Interpreter','latex','FontSize', font_base+2);
    box on
end

function plot_bias_kde(ax, dpp_2d_case, n_shot, clr, font_base, ngrid)
    axes(ax); cla; hold on
    % 1:1 guide (y=1)
    plot([5 55],[1 1],'k-','LineWidth',0.5);
    for i = 1:n_shot
        all3D = 1e9 * dpp_2d_case{i,1,1}; % unfiltered (3D baseline proxy in your earlier usage)
        obs2D = 1e9 * dpp_2d_case{i,3,1}; % spp* = 0.5
        if isempty(all3D) || isempty(obs2D), continue; end
        xmin = 0.9*min([all3D; obs2D]); xmax = 1.1*max([all3D; obs2D]);
        xq = logspace(log10(xmin), log10(xmax), ngrid);
        % Shared bandwidth estimated from the pooled data for stability
        pooled = [all3D; obs2D];
        bw = optimal_kde_bw(pooled);
        f_all = ksdensity(all3D, xq, 'Function','pdf', 'Bandwidth', bw);
        f_obs = ksdensity(obs2D, xq, 'Function','pdf', 'Bandwidth', bw);
        plot(xq, f_obs./max(f_all, eps), 'Color', clr(i,:), 'LineWidth', 2);
    end
    set(gca,'XScale','log','TickLabelInterpreter','latex','FontSize',font_base);
    xlabel('$d_\mathrm{pp}^{(i)}$ [nm]','Interpreter','latex','FontSize',font_base+6);
    xlim([7 50]); ylim([0.6 1.4]);
    title('', 'Interpreter','latex'); % title handled by outer figure context
    box on
end

function bw = optimal_kde_bw(x)
    % Silverman's rule-of-thumb in log domain is more stable for log-spaced x
    x = x(~isnan(x) & isfinite(x) & x>0);
    if numel(x)<2, bw = []; return; end
    lx = log(x);
    sig = std(lx);
    n = numel(lx);
    bw_l = 0.9*min(sig, iqr(lx)/1.34) * n^(-1/5);
    % Map back by approximate scaling: d(log x) ~ dx/x => use multiplicative factor
    % Let bandwidth in x at geometric mean:
    gmx = exp(mean(lx));
    bw = bw_l * gmx; % heuristic; keeps both cases comparable
end

function plot_ratio_vs_da(ax, pars, dpp_2d_case, clr, case_label, font_base)
    axes(ax); cla; hold on
    n_shot = numel(pars);
    % Horizontal guide
    yline(1, 'k-', 'LineWidth', 0.5);
    for i = 1:n_shot
        ratio = dpp_2d_case{i,3,2} ./ dpp_2d_case{i,1,2}; % per-agg GM ratio
        ms_arr = [8 16 12 16 8 8];
        mt_arr = {'^','s','p','*','o','o'};
        scatter(1e9*pars(i).da, ratio, 2*ms_arr(min(i,5)), clr(i,:), mt_arr{min(i,5)}, 'LineWidth', 1);
    end
    set(gca,'XScale','log','TickLabelInterpreter','latex','FontSize',font_base);
    xlabel('$d_\mathrm{a}$ [nm]','Interpreter','latex','FontSize',font_base+6);
    xlim([max(1, 1e9*0.8*min(cat(1,pars.da))), 1e9*1.2*max(cat(1,pars.da))]);
    ylim([0.95 1.16]);
    title(case_label, 'Interpreter','latex','FontSize', font_base+2);
    box on
end

function plot_ratio_vs_snap(ax, pars, dpp_2d_case, clr, xlbl, case_label, font_base, band)
    axes(ax); cla; hold on
    n_shot = numel(pars);
    % Gentle ±band around 1
    if band > 0
        fill([0.25 n_shot+0.75 n_shot+0.75 0.25], [1-band 1-band 1+band 1+band], ...
             [0 0 0], 'FaceAlpha', 0.05, 'EdgeColor','none');
    end
    yline(1, 'k-', 'LineWidth', 0.5);
    for i = 1:n_shot
        r = dpp_2d_case{i,3,2} ./ dpp_2d_case{i,1,2}; % per-agg ratio
        b = boxchart(i*ones(size(r)), r, 'BoxFaceColor', clr(i,:), 'WhiskerLineColor', clr(i,:));
        b.MarkerStyle = '.';
    end
    xticks(1:n_shot); xticklabels(xlbl);
    set(gca,'TickLabelInterpreter','latex','FontSize',font_base);
    xlabel('$n_\mathrm{agg}/(n_\mathrm{agg})_2$ [-]','Interpreter','latex','FontSize',font_base+6);
    ylim([0.95 1.16]);
    title(case_label, 'Interpreter','latex','FontSize', font_base+2);
    box on
end

function plot_sig_ratio_vs_snap(ax, pars, sig2D_case, sig3D_agg, clr, xlbl, case_label, font_base, band)
    axes(ax); cla; hold on
    n_shot = numel(pars);
    if band > 0
        fill([0.25 n_shot+0.75 n_shot+0.75 0.25], [1-band 1-band 1+band 1+band], ...
             [0 0 0], 'FaceAlpha', 0.05, 'EdgeColor','none');
    end
    yline(1, 'k-', 'LineWidth', 0.5);
    for i = 1:n_shot
        r = sig2D_case{i,3} ./ sig3D_agg{i}; % per-agg sigma ratio
        b = boxchart(i*ones(size(r)), r, 'BoxFaceColor', clr(i,:), 'WhiskerLineColor', clr(i,:));
        b.MarkerStyle = '.';
    end
    xticks(1:n_shot); xticklabels(xlbl);
    set(gca,'TickLabelInterpreter','latex','FontSize',font_base);
    xlabel('$n_\mathrm{agg}/(n_\mathrm{agg})_2$ [-]','Interpreter','latex','FontSize',font_base+6);
    ylim([0.95 1.05]);
    title(case_label, 'Interpreter','latex','FontSize', font_base+2);
    box on
end

%% ------------------------ Save outputs & workspace -----------------------

% Output dir
dir0_out = datestr(datetime('now'));
dir0_out = regexprep(dir0_out, ':', '-');
dir0_out = regexprep(dir0_out, ' ', '_');
dir_out = fullfile('outputs', ['PostShield2Case_', dir0_out, filesep]);
if ~isfolder(dir_out), mkdir(dir_out); end

save(fullfile(dir_out, ['PostShield2Case_', dir0_out, '.mat']), ...
     'parsA','parsB','SA','SB','dpp_2d_A','dpp_2d_B','sig_2d_A','sig_2d_B');

exportgraphics(f1, fullfile(dir_out,'Fig1_spp_distributions.jpg'), 'BackgroundColor','none','Resolution',300);
exportgraphics(f2, fullfile(dir_out,'Fig2_sagg_vs_npp.jpg'),      'BackgroundColor','none','Resolution',300);
exportgraphics(f3, fullfile(dir_out,'Fig3_renders_6x3.jpg'),      'BackgroundColor','none','Resolution',300);
exportgraphics(f4, fullfile(dir_out,'Fig4_bias_comp_2x2.jpg'),    'BackgroundColor','none','Resolution',300);
exportgraphics(f5, fullfile(dir_out,'Fig5_ratios_3x2.jpg'),       'BackgroundColor','none','Resolution',300);
