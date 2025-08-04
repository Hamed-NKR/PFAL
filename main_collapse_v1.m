clc
clear
close all

%% initialize

% load previously saved aggregate data

fdir_pars =... % address for folder that contains aggregate coordinate data
    'D:\Hamed\CND\PhD\Publication\DLCA2\mainscatter_sigmapp10\FLAT';
fname_pars = 'FLAT-26NOV24'; % aggregate info filename
parsname = 'pars_out'; % variable to load

fdir_fl =... % address for folder that contains fluid structure data
    'D:\Hamed\CND\PhD\Publication\DLCA2\mainscatter_sigmapp13\SCAT';
fname_fl = 'LD2-25NOV24'; % fluid structure filename
flname = 'fl'; % variable to load

% import library of aggregate information
load(fullfile(fdir_pars, strcat(fname_pars, '.mat')), parsname);

% import fluid data
load(fullfile(fdir_fl, strcat(fname_fl, '.mat')), 'fl');

% change the name of "pars" variable
eval(['pars_in' ' = ' parsname ';']);
eval(['clear ' parsname]);

%% set up effective density figure

% indices of aggregates to be collapsed
ind_agg = [1, 150; 1, 503; 1, 750];

% define colors for plotting
clr1 = hex2rgb({'#DA6C6C', '#659287', '#6096B4'});

% set up effective density figure
f1 = figure(1);
f1.Position = [50, 50, 1200, 700];
set(f1, 'color', 'white')

% create tiled layout to include aggregate renders
tl1 = tiledlayout(3,6);
tl1.TileSpacing = 'compact';
tl1.Padding = 'compact';

tl1_22 = nexttile(1, [3 3]); % make effective density plot tile

% plot constant-mass isobars
dm_cmass = logspace(0, 4, 1000); % define a range of mobility diameters
constMasses = logspace(6, 12, 25); % choose arbitrary constant masses
for k = 1:length(constMasses)
    rho_iso = constMasses(k) ./ (dm_cmass.^3); % proportional to 1/d^3
    plt_cmass = plot(dm_cmass, rho_iso, '--', 'Color', [0.4 0.4 0.4],...
        'LineWidth', 0.5);
    hold on
end

% draw universal correlation line
dm_lim_uc = [1e0 2e4];
n_dm_uc = 1e4;
D_m = 2.48;
rho_eff_100 = 510;
unicor = @(x) rho_eff_100 * (x / 100) .^ (D_m - 3);
r_uc2 = (dm_lim_uc(2) / dm_lim_uc(1)) ^ (1 / (n_dm_uc - 1));
dm_uc = dm_lim_uc(1) * ones(n_dm_uc,1) .* r_uc2 .^ (((1 : n_dm_uc) - 1)');
rho_eff_uc = unicor(dm_uc);
plt_unicor = plot(dm_uc, rho_eff_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);

% draw Sapkota et al.'s universal effective density line (in prep work...
    % ...as of August 2025)
a = 668;
b = 51;
p = 1.7;
colcor = @(x) a + b * (x / 100) .^ (-p);
rho_eff_cc = colcor(dm_uc);
plt_colcor = plot(dm_uc, rho_eff_cc, 'Color', hex2rgb('#113F67'),...
    'LineStyle', ':', 'LineWidth', 2.5);

% set up appearance for effective density plot
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xticks([10:10:90  100:100:900 1000])
xticklabels({'10' '20' '30' '40' '50' '60' '70' '80' '' '100'...
    '200' '300' '400' '500' '600' '700' '800' '' '1000'})
xtickangle(90)
yticks(100:100:1800)
yticklabels({'100' '200' '300' '400' '500' '600' '700' '800' '900' '1000',...
    '1100' '1200' '' '1400' '' '1600' '' '1800'})
xlim([20 800])
ylim([150 1800])
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 24)
ylabel('$\rho_\mathrm{eff} \mathrm{[kg/m^3]}$', 'interpreter', 'latex',...
    'FontSize', 24)

%% simulate collapse and calculate screening

% pick aggregate for collapsing
pp0 = pars_in(ind_agg(1,1)).pp{ind_agg(1,2)};

% apply Beeler et al. (2025)'s collapse algorithm
[pps, n_steps] = PAR.COLAPS(pp0);

% initialize temporal arrays for screening factor (using 3 different...
    % ...calculation methods)

% viewing from one direction
spps_i_singleSide = cell(n_steps, 3);   % per primary particle
spps_singleSide = zeros(n_steps, 3);    % mean

% viewing from both directions
spps_i_doubleSide = cell(n_steps, 3);   % per primary particle
spps_doubleSide = zeros(n_steps, 3);    % mean

% viewing from both directions and neglegting cases where only two...
    % ...primary particles overlap
spps_i_doubleLayer = cell(n_steps, 3);  % per primary particle
spps_doubleLayer = zeros(n_steps, 3);   % mean

% initialize 2d projection of primary particle diameter
dpps_singleSide = zeros(n_steps, 3);
dpps_doubleSide = zeros(n_steps, 3);
dpps_doubleLayer = zeros(n_steps, 3);

% placholders for tracking mass-mobility properties through collapse
dms = zeros(n_steps, 1);      % mobility diameter
rhos = zeros(n_steps, 1);     % effective density

opts_proj.tbar = 'off'; % turn off progress bar for projected area 

% initialize progress bar
disp('Calculate mobility/screening...');
UTILS.TEXTBAR([0, n_steps]);

% calculate effective density, mobility diameter and screening factor...
    % ...during collapse
for kk = 1 : n_steps
    
    % make aggregate strcuture at each timestep (for size calculation)...
        % ...from primary particle data
    pars_temp.pp = pps(kk);
    pars_temp.n = size(pps{kk},1);

    % calculate gyration diameter
    pars_temp = PAR.SIZING(pars_temp);

    % calculate projected area diameter
    pars_temp.da = 2 * sqrt(PAR.PROJECTION(pars_temp, [], 1e2, 5,...
        [], opts_proj)/pi);

    % time history of mobility diameters through collaspe
    dms(kk) = TRANSP.DIAMOBIL(pars_temp.dg(1), pars_temp.da(1), fl);

    % cffective densities as per time history of mobility diameters
    rhos(kk) = 1860 * sum(pps{kk}(:,2).^3) ./ dms(kk).^3;

    % save iteration index where density is near collapse correlation
    if abs(rhos(kk) / colcor(1e9 * dms(kk)) - 1) < 0.01
        n_uni = kk;
    end

    % calculate screening factors for individual primary particles in...
        % ...compacted aggregate
    [spps_i_singleSide{kk}, spps_i_doubleSide{kk},...
        spps_i_doubleLayer{kk}] = SHIELD_METHODS(pps{kk});
    
    % mean screening within each aggregate
    spps_singleSide(kk) = mean(spps_i_singleSide{kk});
    spps_doubleSide(kk) = mean(spps_i_doubleSide{kk});
    spps_doubleLayer(kk) = mean(spps_i_doubleLayer{kk});

    % calculate biased mean primary particle diameter in projection...
        % ...considering screening
    dpps_singleSide(kk) = geomean(pp0(spps_i_singleSide{kk} < 0.5, 2));
    dpps_doubleSide(kk) = geomean(pp0(spps_i_doubleSide{kk} < 0.5, 2));
    dpps_doubleLayer(kk) = geomean(pp0(spps_i_doubleLayer{kk} < 0.5, 2));

    % update progress bar
    UTILS.TEXTBAR([kk, n_steps]);
    
end

%% visualize change in effective density

% indices of intermediate collapsing steps
ind_temp = 2:n_steps-1;
ind_temp(ind_temp == n_uni) = [];

axes(tl1_22) % navigate to effective density subplot

% plot effective densities
scat_intmed = scatter(1e9 * dms(ind_temp), rhos(ind_temp), 5,...
    [0 0 0], '.', 'LineWidth', 1.5); % intermediate
scat_lacey = scatter(1e9 * dms(1), rhos(1), 75, clr1(1,:), 'h',...
    'LineWidth', 2); % initial
scat_unicol = scatter(1e9 * dms(n_uni), rhos(n_uni), 50, clr1(2,:),...
    '^', 'LineWidth', 2); % universal collapse (Sapkota et al.)
pars_final = pars_temp;
scat_colaps = scatter(1e9 * dms(kk), rhos(kk), 50, clr1(3,:),...
    'o', 'LineWidth', 2); % final

% generate legends for various correlations/datasets
legend([plt_unicor, plt_colcor, plt_cmass, scat_intmed, scat_lacey,...
    scat_unicol, scat_colaps], {'Olfert \& Rogak (2019)',...
    'Sapkota et al. (in prep.)', 'Mass isolines', 'Collapsing in process...',...
    'Lacey', 'Mid-collapsed', 'Compact'},...
    'interpreter', 'latex', 'FontSize', 18, 'Location', 'northoutside',...
    'NumColumns', 2)

%% visualize change in aggregate structure

% draw original aggregate
tl1_1_4 = nexttile(4);
opts.cm = UTILS.CUSTOMABYSSMAP('red');
UTILS.PLOTPP(pps{1}(:,3), pps{1}(:,4), pps{1}(:,5), pps{1}(:,2), [], opts)
title({strcat('Aggregate (', int2str(1), ')'), ''},...
    'interpreter', 'latex', 'FontSize', 20)
subtitle('Lacey', 'interpreter', 'latex', 'FontSize', 18)

% draw compacted aggregate to Sapkota et al.'s correlation
tl1_2_4 = nexttile(10);
opts.cm = UTILS.CUSTOMABYSSMAP('green');
UTILS.PLOTPP(pps{n_uni}(:,3), pps{n_uni}(:,4), pps{n_uni}(:,5),...
    pps{n_uni}(:,2), [], opts);
UTILS.SYNC3DVIEW(tl1_1_4, tl1_2_4)
subtitle('Mid-collapsed', 'interpreter', 'latex', 'FontSize', 18)

% draw collapsed aggregate
tl1_3_4 = nexttile(16);
opts.cm = UTILS.CUSTOMABYSSMAP('blue');
UTILS.PLOTPP(pps{kk}(:,3), pps{kk}(:,4), pps{kk}(:,5), pps{kk}(:,2), [],...
    opts);
UTILS.SYNC3DVIEW(tl1_1_4, tl1_3_4)
subtitle('Compact', 'interpreter', 'latex', 'FontSize', 18)

%% visualize change in screening

% set up figure to monitor screening effect
f2 = figure(2);
f2.Position = [100, 100, 900, 300];
set(f2, 'color', 'white')

tl2 = tiledlayout(1,3);
tl2.TileSpacing = 'compact';
tl2.Padding = 'compact';

nexttile(1)

% plot average screening factor vs. mobility diameter
scatter(1e9 * dms(ind_temp), spps_singleSide(ind_temp), 5, [0 0 0], '.',...
    'LineWidth', 1.5);
hold on
scatter(1e9 * dms(1), spps_singleSide(1), 50, clr1(1,:), 'h', 'LineWidth',...
    1.5);
scatter(1e9 * dms(n_uni), spps_singleSide(n_uni), 40, clr1(2,:), '^',...
    'LineWidth', 1.5);
scatter(1e9 * dms(kk), spps_singleSide(kk), 40, clr1(3,:), 'o', 'LineWidth',...
    1.5);
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14,...
    'TickLength', [0.015 0.015])
xlim([0.9*1e9*dms(kk) 1.2*1e9*dms(1)])
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$S_\mathrm{pp}$ [-]', 'interpreter', 'latex',...
    'FontSize', 20)

%% visualize change in primary particle bias

% set up figure to monitor screening effect
f3 = figure(3);
f3.Position = [150, 150, 900, 300];
set(f3, 'color', 'white')

tl3 = tiledlayout(1,3);
tl3.TileSpacing = 'compact';
tl3.Padding = 'compact';

% plot 2d-observed average primary particle diameter vs. mobility diameter
nexttile(1) 
scatter(1e9 * dms(ind_temp), 1e9 * dpps(ind_temp), 5, [0.2 0.2 0.2], '.',...
    'LineWidth', 1.5);
hold on
scatter(1e9 * dms(1), 1e9 * dpps(1), 50, clr1(1,:), 'h', 'LineWidth', 1.5);
scatter(1e9 * dms(n_uni), 1e9 * dpps(n_uni), 40, clr1(2,:), '^',...
    'LineWidth', 1.5);
scatter(1e9 * dms(kk), 1e9 * dpps(kk), 40, clr1(3,:), 'o', 'LineWidth',...
    1.5);
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14,...
    'TickLength', [0.015 0.015])
xlim([0.9*1e9*dms(kk) 1.2*1e9*dms(1)])
ylim([0.9*1e9*min(dpps(1:kk)) 1.2*1e9*max(dpps(1:kk))])
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$d_\mathrm{pp}$ [-]', 'interpreter', 'latex',...
    'FontSize', 20)

%% functionality to calculate shielding factor %%

function [spp_i_singleSide, spp_i_doubleSide, spp_i_doubleLayer] =...
    SHIELD_METHODS(pp, n_points, opts)

% "SHIELD_METHODS": (i) creates equally spaced points on perimeter of...
%   ...each primary particle in z direction, (ii) counts how many other...
%   ...primary particles overlay each of those points on that same...
%   ...direction, and (iii) calculates fraction of overlayed to overall...
%   ...points to decide whether the discretized primary particle is...
%   ...observable or screened based on three different criteria...
%   ...described within the function.
% ----------------------------------------------------------------------- %
%
% Inputs:
%
%   pp: N*6 array of primary particle information for an individual...
%       ...aggregate (column1: index, 2: diameter, 3-5: x,y,z...
%       ...coordinates, 6: subaggregate index).
%
%   n_points: Number of discretization points on each primary particle...
%       ...perimeter (optional, default: 25).
%
%   opts: struct (optional)
%       Options:
%           .chunkSize: Dividing calculations into steps to reduce...
%               ...chances of memory outage (default: 5000).
% ----------------------------------------------------------------------- %
%
% Output:
%
%   spp_i: An N*1 vector where each element is screening factor...
%       ...associated with primary particle of same order in pp.
%
%
%   spp_i in [0 1] range where:
%       spp_i = 1 => screened | spp_i = 0 => observable.
%
%   A discretization point is considered screened if:
%           _singleSide: "One" or more primary particles overlay...
%               ...above it along "+z" projected direction.
%           _doubleSide: "One" or more primary particles overlay with...
%               ...it along both "+/-z" projected directions.
%           _doubleLayer: "Two" or more primary particles overlay...
%               ...with it along both "+/-z" projected directions.
% ----------------------------------------------------------------------- %

% set default values

if nargin < 2 || isempty(n_points); n_points = 25; end

defOpts  = struct('chunkSize', 5000);

if nargin < 3 || isempty(opts)
    opts = defOpts;
else
    fn = fieldnames(defOpts);
    for i = 1:numel(fn)
        if ~isfield(opts, fn{i})
            opts.(fn{i}) = defOpts.(fn{i});
        end
    end
end

% discretization angles
theta = linspace(0, 2*pi, n_points);

npp = size(pp,1); % number of primary particles in aggregate

% initialize screening factor for individual primary particles
spp_i_singleSide = zeros(npp, 1);
spp_i_doubleSide = zeros(npp, 1);
spp_i_doubleLayer = zeros(npp, 1);

% make primary particle pair indices
kk = table2array(combinations(1 : npp, 1 : npp));

% whether a primary particle is further from the viewer than another...
    % ...primary particle
dz = pp(kk(:,1),5) < pp(kk(:,2),5);

% locations of circumferential points on perimeter of each primary...
    % ...particle in x-y plane
x_c = pp(:,3) + repmat(pp(:,2)/2, 1, n_points) .* ...
    cos(repmat(theta, npp, 1)); % x location
y_c = pp(:,4) + repmat(pp(:,2)/2, 1, n_points) .* ...
    sin(repmat(theta, npp, 1)); % y location

for k = 1 : npp % loop through primary particles to calculate...
        % ...screening for each

    % find indices of pair-wise combination specific to each primary...
        % ...particle
    kkk = find((kk(:,1) == k));

    % initialize number of points screened on each primary particle
    spp0 = zeros(1, length(kkk)); 
    
    spp00 = zeros(1, length(kkk)); % considers screening when at least...
        % ... two primary particles fall along the discretization point
    
    n0_points = zeros(1, n_points); % initialize number of primary...
        % ...particles overlaying along each discretization point
    
    % process in chunks to reduce memory usage
    for c = 1 : opts.chunkSize : length(kkk)

        c_end = min(c + opts.chunkSize - 1, length(kkk));
        idx_chunk = kkk(c:c_end);
        src_idx = kk(idx_chunk,1); % always equal to k
        tgt_idx = kk(idx_chunk,2);

        % circumferential point distances to other primary particle...
            % ...centers
        dx = x_c(src_idx,:) - pp(tgt_idx,3);
        dy = y_c(src_idx,:) - pp(tgt_idx,4);
        r2 = dx.^2 + dy.^2;

        % compare to radius of other particle (squared threshold)
        r_thresh2 = pp(tgt_idx,2).^2 / 4;
        r_thresh2 = repmat(r_thresh2, 1, size(dx,2));

        % whether circumferential points fall within projections...
            % ...of other primary particles
        dr_chunk = r2 < r_thresh2;

        % count the number of screened circumferential points
        spp0(c:c_end) = sum(dr_chunk, 2);
        
        % count how many primary particles are along each point
        n0_points = n0_points + sum(dr_chunk, 1);

        % count screened points with the more rigorous layered criteria
        spp00(c:c_end) = sum(dr_chunk(:, n0_points>1), 2);

    end
    
    % select maximum obscured perimeter points as final screening
    if ~isempty(spp0(dz(kkk)))
        spp_i_singleSide(k) = max(spp0(dz(kkk))) / n_points; % choose...
            % ...only the ones with lower z
        spp_i_doubleSide(k) = max(spp0) / n_points; % choose regardless...
            % ...of z
        spp_i_doubleLayer(k) = max(spp00) / n_points; % choose...
            % ...considering a simple optical depth criterion
    end

end

end

