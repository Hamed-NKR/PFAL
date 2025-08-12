%%% Demo: point-touching chain, rotated from z-vertical to horizontal,
%   computing spp under three screening models and plotting a selected
%   orientation with screening factor colorcoding.
%
%   Models (fixed API):
%   - 'opaque'          : front-only (+z), single layer
%   - 'transparent'     : any occluder (ignore z), single layer
%   - 'semitransparent' : any occluder (ignore z), double layer
%
%   Notes:
%   • Perimeter is sampled at n_sample points per primary.

clc; clear; close all;

%% ------------------- USER PARAMETERS ------------------- %%

npp       = 12;          % number of primary particles in the chain
dpp       = 20e-9;       % diameter [m]
n_steps   = 11;          % number of rotation increments (including endpoints)
n_sample  = 64;          % perimeter discretization points for spp

% Rotation path: pitch about +y from 0 to 90 degrees (z --> x)
angles = linspace(0, pi/2, n_steps);

% Which orientation to plot (1..n_steps)
plot_step = 7;           % set to any 1..n_steps

% Screening options
opts_render.chunkSize = 512;   % memory ~ chunkSize * n_sample booleans
opts_render.xyPrune   = true;  % skip targets with center distance > rk + rt

%% ------------------- BUILD BASE CHAIN (ALONG +Z) ------------------- %%

% Centers at z = m*d, centered about origin
idx = (0:npp-1)';
z0  = (idx - (npp-1)/2) * dpp;
x0  = zeros(npp,1);
y0  = zeros(npp,1);

% pp columns convention: [id, diameter, x, y, z, subagg]
pp0 = [(1:npp)', dpp*ones(npp,1), x0, y0, z0, ones(npp,1)];

%% ------------------- STORAGE ------------------- %%

coords = zeros(npp,3,n_steps);         % (x,y,z) per orientation
spp_opaque          = zeros(npp, n_steps);
spp_transparent     = zeros(npp, n_steps);
spp_semitransparent = zeros(npp, n_steps);

%% ------------------- ROTATE, COMPUTE SPP ------------------- %%

for t = 1:n_steps
    a = angles(t);

    % Rotation about +y: R_y(a)
    Ry = [ cos(a)  0  sin(a);
            0      1   0    ;
           -sin(a) 0  cos(a)];

    % Rotate all centers
    XYZ = [pp0(:,3), pp0(:,4), pp0(:,5)] * Ry.';
    pp  = pp0;
    pp(:,3:5) = XYZ;

    % Store coordinates
    coords(:,:,t) = XYZ;

    % Compute spp for all three models
    spp_opaque(:,t)          = SHIELD(pp, n_sample, 'opaque'         , opts_render);
    spp_transparent(:,t)     = SHIELD(pp, n_sample, 'transparent'    , opts_render);
    spp_semitransparent(:,t) = SHIELD(pp, n_sample, 'semitransparent', opts_render);
end

%% ------------------- PLOT SELECTED ORIENTATION ------------------- %%

plot_step = max(1, min(n_steps, round(plot_step)));  % clamp
ang_deg   = angles(plot_step) * 180/pi;

x = coords(:,1,plot_step);
y = coords(:,2,plot_step);
z = coords(:,3,plot_step);
d_vec = pp0(:,2);  % diameter unchanged by rotation

f = figure;
f.Position = [300, 300, 900, 400];
set(f, 'color', 'white')
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% options for rendering
opts_render.cc   = 'on';
opts_render.cm   = colormap('gray'); opts_render.cm = flip(opts_render.cm,1);
opts_render.clim = [0 1];  % enforce full range from 0 to 1
opts_render.ft   = 0.95;

nexttile(1); % opaque
UTILS.PLOTPP_CONTINUOUS(x, y, z, d_vec, spp_opaque(:,plot_step), opts_render);
title(sprintf('opaque (%.1f°)', ang_deg)); axis equal
cb = colorbar('southoutside'); cb.Label.String = 's_{pp}^{(i)} [-]';

nexttile(2); % transparent
UTILS.PLOTPP_CONTINUOUS(x, y, z, d_vec, spp_transparent(:,plot_step), opts_render);
title(sprintf('transparent (%.1f°)', ang_deg)); axis equal
cb = colorbar('southoutside'); cb.Label.String = 's_{pp}^{(i)} [-]';

nexttile(3); % semitransparent
UTILS.PLOTPP_CONTINUOUS(x, y, z, d_vec, spp_semitransparent(:,plot_step), opts_render);
title(sprintf('semitransparent (%.1f°)', ang_deg)); axis equal
cb = colorbar('southoutside'); cb.Label.String = 's_{pp}^{(i)} [-]';

%% =================== LOCAL FUNCTION (STABLE API) ===================

function spp = SHIELD(pp, n_sample, method, opts)
% Compute screening factor for one aggregate under a single method.
%
% Required:
%   pp        : [N x 6] with columns [id, d, x, y, z, subagg]
% Optional (defaults if omitted/empty):
%   n_sample  : perimeter samples per primary (default 64)
%   method    : 'opaque' | 'transparent' | 'semitransparent' (default 'semitransparent')
%   opts      : struct with fields (defaults shown):
%                 .chunkSize (512), .xyPrune (true)
%
% Conventions:
%   Circle radius = d/2. Screening per point:
%     'opaque'          : front-only (+z), threshold >=1 (uses z-tolerance)
%     'transparent'     : ignore z, total threshold >=1
%     'semitransparent' : ignore z, total threshold >=2
%
% Numerics:
%   • strict '<' for non-concentric pairs (tangency doesn’t count)
%   • concentric occluders (co-centered within tol) cover the *entire* perimeter
%   • front/back test uses a small z-tolerance

    % ---- defaults ----
    if nargin < 2 || isempty(n_sample), n_sample = 64; end
    if nargin < 3 || isempty(method),    method    = 'semitransparent'; end
    def.chunkSize = 512; def.xyPrune = true;
    if nargin < 4 || isempty(opts), opts = def; else
        if ~isfield(opts,'chunkSize') || isempty(opts.chunkSize), opts.chunkSize = def.chunkSize; end
        if ~isfield(opts,'xyPrune')   || isempty(opts.xyPrune),   opts.xyPrune   = def.xyPrune;   end
    end

    npp = size(pp,1);
    spp = zeros(npp,1,'double');

    d  = pp(:,2);  r  = 0.5*d;
    xc = pp(:,3);  yc = pp(:,4);  zc = pp(:,5);

    % Perimeter sampling (uniform, no half-step)
    ang = (0:n_sample-1) * (2*pi / max(1,n_sample));

    % Tolerances
    tolz   = 1e-12 * max(1, max(abs(zc)) + max(r));  % z-order tol
    tolc   = 1e-12 * max(1, max(r));                 % co-center tol (xy)

    for k = 1:npp
        rk = r(k);
        xk = xc(k) + rk * cos(ang);
        yk = yc(k) + rk * sin(ang);

        if npp <= 1, spp(k) = 0; continue; end
        tgt_all = setdiff(1:npp, k);

        % Optional safe XY pruning
        if opts.xyPrune && ~isempty(tgt_all)
            dx0  = xc(tgt_all) - xc(k);
            dy0  = yc(tgt_all) - yc(k);
            rtgt = r(tgt_all);
            keep = (dx0.^2 + dy0.^2) <= (rk + rtgt).^2;
            tgt_all = tgt_all(keep);
            if isempty(tgt_all), spp(k) = 0; continue; end
        end

        % Initialize per-point counters
        switch method
            case 'opaque'
                n_plus = zeros(1, n_sample, 'uint16');   % +z occluders
            case 'transparent'
                n_all  = zeros(1, n_sample, 'uint16');   % total (ignore z)
            case 'semitransparent'
                n_all  = zeros(1, n_sample, 'uint16');   % total (ignore z)
            otherwise
                error('Unknown method: %s', method);
        end

        % Stream targets in chunks (memory-safe)
        cs = min(opts.chunkSize, max(1, numel(tgt_all)));
        for c = 1:cs:numel(tgt_all)
            idx = tgt_all(c : min(c+cs-1, numel(tgt_all)));

            xt = xc(idx);
            yt = yc(idx);
            zt = zc(idx);
            rt = r(idx);

            % 1) Handle concentric occluders first (cover entire perimeter)
            dist_xy = hypot(xt - xc(k), yt - yc(k));
            conc    = dist_xy <= tolc;  % co-centered (within tol)

            if any(conc)
                switch method
                    case 'opaque'
                        front_conc = conc & ((zt - zc(k)) > tolz);
                        if any(front_conc)
                            n_plus = n_plus + uint16(sum(front_conc)) * ones(1, n_sample, 'uint16');
                        end
                    case 'transparent'
                        n_all = n_all + uint16(sum(conc)) * ones(1, n_sample, 'uint16');
                    case 'semitransparent'
                        n_all = n_all + uint16(sum(conc)) * ones(1, n_sample, 'uint16');
                end
            end

            % 2) Non-concentric occluders: strict '<' so tangency doesn’t count
            non = ~conc;
            if any(non)
                xt = xt(non); yt = yt(non); zt2 = zt(non); rt = rt(non);
                dx = xk - xt;                 % [numel(non) x n_sample]
                dy = yk - yt;
                inside = (dx.^2 + dy.^2) < (rt.^2);  % strict

                switch method
                    case 'opaque'
                        front = (zt2 - zc(k)) > tolz;
                        if any(front)
                            n_plus = n_plus + uint16(sum(inside(front,:), 1));
                        end
                    case 'transparent'
                        n_all = n_all + uint16(sum(inside, 1));
                    case 'semitransparent'
                        n_all = n_all + uint16(sum(inside, 1));
                end
            end
        end

        switch method
            case 'opaque'
                spp(k) = mean(n_plus >= 1);
            case 'transparent'
                spp(k) = mean(n_all  >= 1);
            case 'semitransparent'
                spp(k) = mean(n_all  >= 2);
        end
    end
end
