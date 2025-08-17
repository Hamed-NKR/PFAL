%%% Demo: point-touching chain, rotated from z-vertical to horizontal,
%   computing spp under three screening models and plotting a selected
%   orientation with screening factor colorcoding.
%
%   Models:
%   - 'opaque': front-only (+z), single layer
%   - 'transparent': any occluder (ignore z), single layer
%   - 'semitransparent': any occluder (ignore z), double layer
%
%   Notes:
%   - Perimeter is sampled at n_sample points per primary.

clc; clear; close all;

%% USER PARAMETERS

npp       = 10;          % number of primary particles in the chain
dpp       = 20e-9;       % diameter [m]
n_steps   = 4;          % number of rotation increments (including endpoints)
n_sample  = 64;          % perimeter discretization points for spp

% Rotation path: pitch about +y from 0 to 90 degrees (z --> x)
angles = linspace(0, pi/2, n_steps);

% Screening options
opts_render.chunkSize = 512;   % memory ~ chunkSize * n_sample booleans
opts_render.xyPrune   = true;  % skip targets with center distance > rk + rt

%% BUILD BASE CHAIN (ALONG +Z)

% Centers at z = m*d, centered about origin
idx = (0:npp-1)';
z0  = (idx - (npp-1)/2) * dpp;
x0  = zeros(npp,1);
y0  = zeros(npp,1);

% pp columns convention: [id, diameter, x, y, z, subagg]
pp0 = [(1:npp)', dpp*ones(npp,1), x0, y0, z0, ones(npp,1)];

%% STORAGE

coords = zeros(npp,3,n_steps);         % (x,y,z) per orientation
spp_opaque          = zeros(npp, n_steps);
spp_transparent     = zeros(npp, n_steps);
spp_semitransparent = zeros(npp, n_steps);

%% ROTATE, COMPUTE SPP

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
    spp_opaque(:,t)          = PAR.SHIELD(pp, n_sample, 'opaque'         , opts_render);
    spp_transparent(:,t)     = PAR.SHIELD(pp, n_sample, 'transparent'    , opts_render);
    spp_semitransparent(:,t) = PAR.SHIELD(pp, n_sample, 'semitransparent', opts_render);
end

%% PLOT SELECTED ORIENTATION

f = figure;
f.Position = [50, 50, 600, 900];
set(f, 'color', 'white')

% Which orientations to plot
if n_steps >= 4
    % generate 4 near-equally spaced orientations between 1 and n_steps
    ind_plot = round(linspace(1, n_steps, 4));
else
    % single middle-range orientation
    ind_plot = round((n_steps + 1) / 2);
end

% assign how many rows exist in tiled layout
n_plot = length(ind_plot);
tiledlayout(n_plot,3,'TileSpacing','compact','Padding','compact');

% options for rendering
opts_render.cc   = 'on';
opts_render.cm   = colormap('gray'); opts_render.cm = flip(opts_render.cm,1);
opts_render.clim = [0 1];  % enforce full range from 0 to 1
opts_render.ft   = 1;

for i = 1 : n_plot
    
    % find primary particle coordinates for each orientations
    ang_deg   = angles(ind_plot(i)) * 180/pi;
    x = coords(:,1,ind_plot(i));
    y = coords(:,2,ind_plot(i));
    z = coords(:,3,ind_plot(i));
    d_vec = pp0(:,2); % diameter unchanged by rotation

    nexttile % opaque
    UTILS.PLOTPP_CONTINUOUS(x, y, z, d_vec, spp_opaque(:,ind_plot(i)),...
        opts_render);
    axis equal
    if i == 1
        title('Opaque', 'interpreter', 'latex', 'FontSize', 18);
    end
    subtitle(' ', 'interpreter', 'latex', 'FontSize', 16);

    nexttile % transparent
    UTILS.PLOTPP_CONTINUOUS(x, y, z, d_vec, spp_transparent(:,ind_plot(i)),...
        opts_render);
    axis equal
    if i == 1
        title('Transparent', 'interpreter', 'latex', 'FontSize', 18);
    end
    subtitle(sprintf('($%.1f^{\\circ}$)', ang_deg), 'interpreter',...
        'latex', 'FontSize', 16);
    
    % draw colormap
    if i == n_plot
        cb = colorbar('southoutside');
        cb.FontSize = 14;
        cb.Label.FontSize = 20;
        cb.TickLabelInterpreter = 'latex';
        cb.Label.Interpreter = 'latex';
        cb.Label.String = '$S_\mathrm{pp}^\mathrm{(i)}$ [-]';
        cb.LineWidth = 1;
        cb.TickLength = 0.02;
    end
    
    nexttile % semitransparent
    UTILS.PLOTPP_CONTINUOUS(x, y, z, d_vec,...
        spp_semitransparent(:,ind_plot(i)), opts_render);
    axis equal
    if i == 1
        title('Semitransparent', 'interpreter', 'latex', 'FontSize', 18);
    end
    subtitle(' ', 'interpreter',...
        'latex', 'FontSize', 16);

end