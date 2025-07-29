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

% set up effective density figure
f1 = figure(1);
f1.Position = [50, 50, 900, 600];
set(f1, 'color', 'white')

% create tiled layout to include aggregate renders
tl1 = tiledlayout(2,3);
tl1.TileSpacing = 'compact';
tl1.Padding = 'compact';

tl1_22 = nexttile(1, [2 2]); % make effective density plot tile

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

% draw Nishan Sapkota's universal effective density line (in prep work...
    % as of July 2025)
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
    '1100' '1200' '' '1400' '' '' '' '1800'})
xlim([20 800])
ylim([150 1800])
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 24)
ylabel('$\rho_\mathrm{eff} \mathrm{[kg/m^3]}$', 'interpreter', 'latex',...
    'FontSize', 24)

% pick one aggregate
pp0 = pars_out(1).pp{503};
d0 = pp0(:,2); % diameters
x0 = [pp0(:,3), pp0(:,4), pp0(:,5)]; % positions

% draw original aggregate
tl1_1_3 = nexttile(3);
opts.cm = UTILS.CUSTOMABYSSMAP('orange');
UTILS.PLOTPP(x0(:,1), x0(:,2), x0(:,3), d0, [], opts)

% normalize radius and center position
xcm0 = sum((d0 .^ 3) .* x0, 1) ./ sum(d0 .^ 3);
r0 = d0 / 2;
rm0 = mean(r0);
r = r0 / rm0;
x = x0 / rm0;
x = x - mean(x,1);

dtemp = d0;

% collapse parameters
k = 0.2;        % spring to center
k_decay = 1.0;  % how spring weakens
eps = 0.5;      % well depth of LJ
v_decay = 0.98; % friction
jit = 0.0;      % random jitter
jit_decay = 0.98;
Dt = 0.01;      % time step
nt = 1000;      % steps per block
nj = 10;         % number of blocks

% precompute values
N = size(x,1);
sig = (r + r') / (2^(1/6)); % equilibrium distance
v = zeros(N,3);             % initial velocities
m = (4/3) * pi * r.^3;      % masses
ii=1;

% initialize temporal array for screening factor
spp_i = cell(nj*nt+1, 1); % per primary particle
spp = zeros(nj*nt+1, 1); % mean

% calculate screening in each primary particle at initial moment
spp_i{1} = SHIELD(pp0);

spp(1) = mean(spp_i{1}); % calculate average per aggregate at time 0

% initialize 2d projection of primary particle diameter
dpp = zeros(nj*nt+1, 1);
dpp(1) = geomean(pp0(spp_i{1} < 0.5, 2)); % average between primary...
    % ...particles with screening below 0.5

dm = zeros(nj*nt+1, 1); % placeholder for intermediate mobility diameters
rho_eff = zeros(nj*nt+1, 1); % placeholder for intermediate effective densities

stopit = false; % flag for terminating iterations

opts_proj.tbar = 'off'; % turn off progress bar for projected area 

% initialize progress bar
disp('Aggregate collapsing...');
UTILS.TEXTBAR([0, nj * nt]);

lastInd = 1; % index for last iteration number

% simulate collapse (adopting the simplified algorithm proposed by...
    % ...Beeler et al. (2025): "A Framework for Quantifying the Size and...
    % ...Fractal Dimension of Compacting Soot Particles")
for jj = 1:nj % repeat the entire integration sequence (number of blocks)

    for ii = 1:nt % number of time integration steps in each block...
        % ...(e.g. for coating)

        % apply decay to jitter and spring for numerical stability
        jit = jit * jit_decay;
        k = k * k_decay;

        % center of mass force
        cen = sum(x .* r, 1) ./ sum(r);
        F = -k * (x - cen);

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
        Fvw_mag = 48 * eps * ((sig ./ d0_clamped).^12 ./ d0_clamped - ...
                              0.5 * (sig ./ d0_clamped).^6 ./ d0_clamped);
        Fvw_mag(dtemp > 1.5 .* sig) = 0;
        
        % apply force direction
        Fvw_mag_exp = repmat(Fvw_mag, 1, 1, 3);
        Fvw = squeeze(sum(Fvw_mag_exp .* d_unit, 2));
        F = F + Fvw + jit * randn(N,3);

        % leapfrog integration
        a = F ./ m;
        v = v + a * Dt;
        x = x + v * Dt;
        v = v * v_decay;

        % scale aggregate back to original coordinates for strcutural...
            % ...analysis
        x_temp = rm0 * x;
        xcm_temp = sum((d0 .^ 3) .* x_temp, 1) ./ sum(d0 .^ 3);
        x_temp = x_temp + (xcm0 - xcm_temp);

        % calculate temporal effective density
        pp_temp = pp0;
        pp_temp(:,3:5) = x_temp;
        pars_temp.pp = {pp_temp};
        pars_temp.n = size(pp_temp,1);
        pars_temp = PAR.SIZING(pars_temp);
        pars_temp.da = 2 * sqrt(PAR.PROJECTION(pars_temp, [], 1e2, 5,...
            [], opts_proj)/pi);
        dm(lastInd) = TRANSP.DIAMOBIL(pars_temp.dg(1), pars_temp.da(1), fl);
        rho_eff(lastInd) = 1860 * sum(pars_temp.pp{1}(:,2).^3) ./...
            dm(lastInd).^3;
        
        % update progress bar
        UTILS.TEXTBAR([ii + (nt*(jj-1)), nj*nt]);
        
        % calculate screening for compacted aggregate
        spp_i{lastInd} = SHIELD(pp_temp);
        spp(lastInd) = mean(spp_i{lastInd});
        
        % calculate biased projected mean primary particle diameter
        dpp(lastInd) = geomean(pp0(spp_i{lastInd} < 0.5, 2));
        
        % decide on terminating iterations based on criteria on...
            % ...effective density
        if rho_eff(lastInd) > 10 * colcor(1e9 * dm(lastInd)) ||...
                lastInd == nj * nt
            stopit = true; % for outer loop (used later)
            break   % exits the inner loop
            
        else
            lastInd = lastInd + 1; % update iteration index
        end

    end

    if stopit
        break   % exits the outer loop
    end

end

axes(tl1_22) % navigate to effective density subplot for intermediate...
    % ...collapse process data monitoring

% mark intermediate effective density
scat_intmed = scatter(1e9 * dm(2:lastInd-1), rho_eff(2:lastInd-1),...
    5, [0 0 0], '.', 'LineWidth', 1.5);

% calculate and mark initial effective density
dm(1) = TRANSP.DIAMOBIL(pars_out(1).dg(503), pars_out(1).da(503), fl);
rho_eff(1) = 1860 * sum(pars_out(1).pp{503}(:,2).^3) ./ dm(1).^3;
scat_lacey = scatter(1e9 * dm(1), rho_eff(1), 75, hex2rgb('#E16A54'), 'h',...
    'LineWidth', 2);

% mark final effective density
pars_final = pars_temp;
scat_colaps = scatter(1e9 * dm(lastInd), rho_eff(lastInd), 40,...
    hex2rgb('#0ABAB5'), 'o', 'LineWidth', 2);

% generate legends for effective density plot
legend([plt_unicor, plt_colcor, plt_cmass, scat_lacey, scat_intmed,...
    scat_colaps], {'Olfert \& Rogak (2019)', 'Sapkota et al. (in prep.)',...
    'Mass isolines', 'Lacey aggregate', 'Collapsing...', 'Compact aggregate'},...
    'interpreter', 'latex', 'FontSize', 14, 'Location', 'southwest')

% draw collapsed aggregate
tl1_2_3 = nexttile(6);
opts.cm = UTILS.CUSTOMABYSSMAP('blue');
UTILS.PLOTPP(x_temp(:,1), x_temp(:,2), x_temp(:,3), d0, [], opts);
UTILS.SYNC3DVIEW(tl1_1_3, tl1_2_3)

% set up figure to monitor screening effect
f2 = figure(2);
f2.Position = [100, 100, 800, 400];
set(f2, 'color', 'white')

tl2 = tiledlayout(1,2);
tl2.TileSpacing = 'compact';
tl2.Padding = 'compact';

nexttile(1) % plot average screening factor vs. mobility diameter
scatter(1e9 * dm(2:lastInd-1), spp(2:lastInd-1), 5, [0 0 0], '.',...
    'LineWidth', 1.5);
hold on
scatter(1e9 * dm(1), spp(1), 45, hex2rgb('#E16A54'), 'h',...
    'LineWidth', 1.5);
scatter(1e9 * dm(lastInd), spp(lastInd), 30, hex2rgb('#0ABAB5'), 'o',...
    'LineWidth', 1.5);
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14,...
    'TickLength', [0.015 0.015])
xlim([0.9*1e9*dm(lastInd) 1.2*1e9*dm(1)])
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$S_\mathrm{pp}$ [-]', 'interpreter', 'latex',...
    'FontSize', 20)

nexttile(2) % plot 2d-observed average primary particle diameter...
    % ...vs. mobility diameter
scatter(1e9 * dm(2:lastInd-1), 1e9 * dpp(2:lastInd-1), 5, [0 0 0], '.',...
    'LineWidth', 1.5);
hold on
scatter(1e9 * dm(1), 1e9 * dpp(1), 45, hex2rgb('#E16A54'), 'h',...
    'LineWidth', 1.5);
scatter(1e9 * dm(lastInd), 1e9 * dpp(lastInd), 30, hex2rgb('#0ABAB5'),...
    'o', 'LineWidth', 1.5);
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14,...
    'TickLength', [0.015 0.015])
xlim([0.9*1e9*dm(lastInd) 1.2*1e9*dm(1)])
ylim([0.9*1e9*min(dpp(1:lastInd)) 1.2*1e9*max(dpp(1:lastInd))])
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$d_\mathrm{pp}$ [-]', 'interpreter', 'latex',...
    'FontSize', 20)

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