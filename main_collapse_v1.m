clc
clear
close all

% load previously saved aggregate data
fdir_in =... % address for folder that contains aggregate coordinate data
    'D:\Hamed\CND\PhD\Publication\DLCA2\mainscatter_sigmapp10\FLAT';
fname_in = 'FLAT-26NOV24'; % aggregate info filename
varnames = {'pars_out'}; % variable to load
for i = 1 : numel(varnames)
    % import library of aggregate information
    load(fullfile(fdir_in, strcat(fname_in, '.mat')), varnames{i});

    % import fluid data
    load(fullfile('D:\Hamed\CND\PhD\Publication\DLCA2\mainscatter_sigmapp13\SCAT',...
        strcat('LD2-25NOV24', '.mat')), 'fl');    
end

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
    % text(dm_cmass(end)*0.95, rho_iso(end), ...
    %     sprintf('M = %.0e', constMasses(k)), ...
    %     'FontSize', 9, 'HorizontalAlignment', 'right');
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

% pick aggregate for collapsing
pp0 = pars_out(1).pp{503};

% apply Beeler et al. (2025)'s collapse algorithm
[pps, n_steps] = COLAPS(pp0);

% initialize temporal array for screening factor
spps_i = cell(n_steps, 1); % per primary particle
spps = zeros(n_steps, 1); % mean

% calculate screening in each primary particle at initial moment
spps_i{1} = SHIELD(pp0);

spps(1) = mean(spps_i{1}); % calculate average per aggregate at time 0

% initialize 2d projection of primary particle diameter
dpps = zeros(n_steps, 1);
dpps(1) = geomean(pp0(spps_i{1} < 0.5, 2)); 

opts_proj.tbar = 'off'; % turn off progress bar for projected area 

% placholders for tracking mass-mobility properties through collapse
dms = zeros(n_steps, 1);      % mobility diameter
rhos = zeros(n_steps, 1);     % effective density

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
    spps_i{kk} = SHIELD(pps{kk});
    spps(kk) = mean(spps_i{kk}); % mean screening within each aggregate

    % calculate biased mean primary particle diameter in projection...
        % ...considering screening
    dpps(kk) = geomean(pp0(spps_i{kk} < 0.5, 2));

end

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
    'interpreter', 'latex', 'FontSize', 16, 'Location', 'northoutside',...
    'NumColumns', 2)

% draw original aggregate
tl1_1_4 = nexttile(4);
opts.cm = UTILS.CUSTOMABYSSMAP('red');
UTILS.PLOTPP(pps{1}(:,3), pps{1}(:,4), pps{1}(:,5), pps{1}(:,2), [], opts)
title({cell2mat(strcat('Aggregate', {' '}, int2str(1))), ''},...
    'interpreter', 'latex', 'FontSize', 20)
subtitle('Lacey', 'interpreter', 'latex', 'FontSize', 16)


% draw compacted aggregate to Sapkota et al.'s correlation
tl1_2_4 = nexttile(10);
opts.cm = UTILS.CUSTOMABYSSMAP('green');
UTILS.PLOTPP(pps{n_uni}(:,3), pps{n_uni}(:,4), pps{n_uni}(:,5),...
    pps{n_uni}(:,2), [], opts);
UTILS.SYNC3DVIEW(tl1_1_4, tl1_2_4)
subtitle('Mid-collapsed', 'interpreter', 'latex', 'FontSize', 16)

% draw collapsed aggregate
tl1_3_4 = nexttile(16);
opts.cm = UTILS.CUSTOMABYSSMAP('blue');
UTILS.PLOTPP(pps{kk}(:,3), pps{kk}(:,4), pps{kk}(:,5), pps{kk}(:,2), [], opts);
UTILS.SYNC3DVIEW(tl1_1_4, tl1_3_4)
subtitle('Compact', 'interpreter', 'latex', 'FontSize', 16)

% set up figure to monitor screening effect
f2 = figure(2);
f2.Position = [100, 100, 800, 400];
set(f2, 'color', 'white')

tl2 = tiledlayout(1,2);
tl2.TileSpacing = 'compact';
tl2.Padding = 'compact';

nexttile(1) % plot average screening factor vs. mobility diameter
scatter(1e9 * dms(ind_temp), spps(ind_temp), 5, [0 0 0], '.',...
    'LineWidth', 1.5);
hold on
scatter(1e9 * dms(1), spps(1), 50, clr1(1,:), 'h', 'LineWidth', 1.5);
scatter(1e9 * dms(n_uni), spps(n_uni), 40, clr1(2,:), '^', 'LineWidth', 1.5);
scatter(1e9 * dms(kk), spps(kk), 40, clr1(3,:), 'o', 'LineWidth', 1.5);
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14,...
    'TickLength', [0.015 0.015])
xlim([0.9*1e9*dms(kk) 1.2*1e9*dms(1)])
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$S_\mathrm{pp}$ [-]', 'interpreter', 'latex',...
    'FontSize', 20)

nexttile(2) % plot 2d-observed average primary particle diameter...
    % ...vs. mobility diameter
scatter(1e9 * dms(ind_temp), 1e9 * dpps(ind_temp), 5, [0 0 0], '.',...
    'LineWidth', 1.5);
hold on
scatter(1e9 * dms(1), 1e9 * dpps(1), 50, clr1(1,:), 'h', 'LineWidth', 1.5);
scatter(1e9 * dms(n_uni), 1e9 * dpps(n_uni), 40, clr1(2,:), '^',...
    'LineWidth', 1.5);
scatter(1e9 * dms(kk), 1e9 * dpps(kk), 40, clr1(3,:), 'o', 'LineWidth', 1.5);
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14,...
    'TickLength', [0.015 0.015])
xlim([0.9*1e9*dms(kk) 1.2*1e9*dms(1)])
ylim([0.9*1e9*min(dpps(1:kk)) 1.2*1e9*max(dpps(1:kk))])
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$d_\mathrm{pp}$ [-]', 'interpreter', 'latex',...
    'FontSize', 20)

%% Function to collapse %%

function [pps, n_steps] = COLAPS(pp0, n_steps, specs, opts)
% "COLAPS" compacts DLCA aggregates using spring and van der Waals...
%   ...forces based on Beeler et al. (2025): "A Framework for...
%   ...Quantifying the Size and Fractal Dimension of Compacting Soot...
%   ...Particles".
% ----------------------------------------------------------------------- %
%
% Inputs:
%
%   pp0: N*6 array of primary particle information for an individual...
%       ...aggregate (column1: index, 2: diameter, 3-5: x,y,z...
%       ...coordinates, 6: subaggregate index).
%
%   n_steps: Number of iterations (optional, default: 10,000)
%
%   specs: struct (optional)
%       Simulation parameters:
%           .k_spring   : Spring constant (default: 0.2)
%           .k_decay    : Spring weakening rate (default: 1.0)
%           .jit        : Random perturbation factor (default: 0.01)
%           .jit_decay  : perturbation weakening rate (default: 0.98)
%           .lj_eps     : Lennard-Jones well depth (default: 0.5)
%           .lj_decay   : Lennard-Jones sigma (default: 0.98)
%           .dt         : Timestep (default: 0.01)
%
%   opts: struct (optional)
%       Options:
%           .showProgress : Display progress bar (default: true)
% ----------------------------------------------------------------------- %
%
% Output:
%
%   pps: A cell array containing time history of aggregate's primary...
%       ...particle information over the course of collaspe.
% ----------------------------------------------------------------------- %

% default values
defSpecs = struct('k_spring', 0.2, 'k_decay', 1.0, 'jit', 0.01,...
    'jit_decay', 0.98, 'lj_eps', 0.5, 'lj_decay', 0.98, 'dt', 0.01);
defOpts  = struct('showProgress', true);

% merge user inputs with defaults

if nargin < 2 || isempty(n_steps); n_steps = 1e4; end

if nargin < 3 || isempty(specs)
    specs = defSpecs;
else
    fn = fieldnames(defSpecs);
    for i = 1:numel(fn)
        if ~isfield(specs, fn{i})
            specs.(fn{i}) = defSpecs.(fn{i});
        end
    end
end

if nargin < 4 || isempty(opts)
    opts = defOpts;
else
    fn = fieldnames(defOpts);
    for i = 1:numel(fn)
        if ~isfield(opts, fn{i})
            opts.(fn{i}) = defOpts.(fn{i});
        end
    end
end

d0 = pp0(:,2); % extract primary particle diameters
x0 = [pp0(:,3), pp0(:,4), pp0(:,5)]; % primary particle [x y z] coordinates

% normalize radius and center position
xcm0 = sum((d0 .^ 3) .* x0, 1) ./ sum(d0 .^ 3);
r0 = d0 / 2;
rm0 = mean(r0);
r = r0 / rm0;
x = x0 / rm0;
x = x - mean(x,1);

pps = cell(n_steps, 1);    % initialize primary particle coordinates
npp = size(x,1);           % number of primary particles

% precompute values
sig = (r + r') / (2^(1/6));   % equilibrium distance
v = zeros(npp,3);             % initial velocities
m = (4/3) * pi * r.^3;        % masses

if opts.showProgress
    % initialize progress bar for collapse
    disp('Aggregate collapsing...');
    UTILS.TEXTBAR([0, n_steps]);
end

kk = 1; % index for iteration

% simulate collapse, adopting Beeler et al. (2025)'s algorithm
while kk <= n_steps
        
        % apply decay to jitter and spring for numerical stability
        specs.jit = specs.jit * specs.jit_decay;
        specs.k_spring = specs.k_spring * specs.k_decay;

        % center of mass force
        cen = sum(x .* r, 1) ./ sum(r);
        F = -specs.k_spring * (x - cen);

        % pairwise distances and unit vectors
        dtemp = squareform(pdist(x));
        dtemp(dtemp == 0) = Inf;
        dvec = permute(x, [1 3 2]) - permute(x, [3 1 2]);
        ds = sqrt(sum(dvec.^2, 3));
        ds(ds == 0) = Inf;
        d_unit = dvec ./ ds;

        % van der waals force
        rmin = 0.3;
        d0_clamped = max(dtemp, rmin);
        Fvw_mag = 48 * specs.lj_eps * ((sig ./ d0_clamped).^12 ./...
            d0_clamped - 0.5 * (sig ./ d0_clamped).^6 ./...
            d0_clamped);
        Fvw_mag(dtemp > 1.5 .* sig) = 0;
        
        % apply force direction
        Fvw_mag_exp = repmat(Fvw_mag, 1, 1, 3);
        Fvw = squeeze(sum(Fvw_mag_exp .* d_unit, 2));
        F = F + Fvw + specs.jit * randn(npp,3);

        % leapfrog integration
        a = F ./ m;
        v = v + a * specs.dt;
        x = x + v * specs.dt;
        v = v * specs.lj_decay;

        % scale aggregate back to original coordinates for strcutural...
            % ...analysis
        x_temp = rm0 * x;
        xcm_temp = sum((d0 .^ 3) .* x_temp, 1) ./ sum(d0 .^ 3);
        x_temp = x_temp + (xcm0 - xcm_temp);

        % calculate temporal effective density
        pps{kk} = pp0;
        pps{kk}(:,3:5) = x_temp;
        
        if opts.showProgress
            UTILS.TEXTBAR([kk, n_steps]); % update progress bar
        end
        
        kk = kk + 1; % update iteration index        

end 
    
end

%% Function to calculate shielding factor %%
function spp_i = SHIELD(pp)

n_shield = 25; % resolution for perimeter discretization

% discretize orientation around primary particles for shielding calculation
theta = linspace(0, 2*pi, n_shield);

npp = size(pp,1); % number of primary particles in aggregate

% initialize shielding factor for individual primary particles
spp_i = zeros(npp, 1);

% calculate shielding factor (with spp: [0 1] where spp = 1...
    % ...meaning fully screened) for each primary particle in...
    % ...each aggregate

% make primary particle pair indices
kk = table2array(combinations(1 : npp, 1 : npp));

% whether a primary particle is further from the viewer than...
    % ...another primary particle
dz = pp(kk(:,1),5) < pp(kk(:,2),5);

% locations of circumferential points on perimeter of each...
% ...primary particle in x-y plane

% x location of circumferential points
x_c = pp(:,3) + repmat(pp(:,2)/2, 1, n_shield) .* ...
    cos(repmat(theta, npp, 1));

% y location
y_c = pp(:,4) + repmat(pp(:,2)/2, 1, n_shield) .* ...
    sin(repmat(theta, npp, 1));

% set the maximum amount of shielding for each primary particle...
% ...from all other primary particles as the final value...
% ...assigned for shielding

chunkSize = 5000; % adjust this based on available memory

for k = 1 : npp

    % find indices of pair-wise combination specific to each...
    % ...primary particle
    kkk = find((kk(:,1) == k) & dz); % choose only the ones...
    % ...with lower z

    % delete raw values of number of primary particles shielded...
    % ...from previous iteration
    spp0 = zeros(1, length(kkk));

    % process in chunks to reduce memory usage
    for c = 1 : chunkSize : length(kkk)
        c_end = min(c + chunkSize - 1, length(kkk));
        idx_chunk = kkk(c:c_end);
        src_idx = kk(idx_chunk,1); % always equal to k
        tgt_idx = kk(idx_chunk,2);

        % circumferential point distances to other primary...
        % ...particle centers
        dx = x_c(src_idx,:) - pp(tgt_idx,3);
        dy = y_c(src_idx,:) - pp(tgt_idx,4);
        r2 = dx.^2 + dy.^2;

        % compare to radius of other particle (squared threshold)
        r_thresh2 = pp(tgt_idx,2).^2 / 4;
        r_thresh2 = repmat(r_thresh2, 1, size(dx,2));

        % whether circumferential points fall within projections...
        % ...of other primary particles
        dr_chunk = r2 < r_thresh2;

        % count the number of shielded circumferential points
        spp0(c:c_end) = sum(dr_chunk, 2);
    end

    % select maximum value as final shielding
    if ~isempty(spp0)
        spp_i(k) = max(spp0) / n_shield;
    end

end

end