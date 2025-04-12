function h_pp = PLOTPP_CONTINUOUS(x_pp, y_pp, z_pp, d_pp, colorval, opts)
% PLOTPP_CONTINUOUS displays primary particles in 3D space with scalar or
% categorical coloring.
%
% Inputs:
%   x_pp, y_pp, z_pp: Coordinates of primary particles
%   d_pp: Diameters of primary particles
%   colorval: Scalar values (e.g., shielding) or group IDs for coloring
%   opts: Struct with fields:
%       - cm: colormap (default = gray)
%       - ft: face transparency (default = 1)
%       - cc: 'on' for scalar coloring, 'off' for group-based coloring
%       - cloc: color center location (used for group-based coloring)
%       - clim: [min max] color limits for scalar coloring (optional)
%
% Output:
%   h_pp: Handles to surface objects for each primary particle

% Number of particles
n_pp = numel(x_pp);

% Default options
if ~exist('opts', 'var'), opts = struct(); end
if ~isfield(opts, 'cm'), opts.cm = gray; end
if ~isfield(opts, 'ft'), opts.ft = 1; end
if ~isfield(opts, 'cc'), opts.cc = 'off'; end
if ~isfield(opts, 'cloc'), opts.cloc = 0.25; end

% Determine coloring mode
use_scalar_coloring = strcmpi(opts.cc, 'on') && isnumeric(colorval) && ...
                      (any(colorval ~= round(colorval)) || any(colorval < 1));

% Prepare color mapping
if use_scalar_coloring
    % Use user-defined clim or compute from data
    if isfield(opts, 'clim') && numel(opts.clim) == 2
        vmin = opts.clim(1);
        vmax = opts.clim(2);
    else
        vmin = min(colorval);
        vmax = max(colorval);
    end

    colormap(opts.cm);       % assign colormap to axes
    clim([vmin vmax]);      % set color limits
else
    ids = unique(colorval);
    n_ids = length(ids);
    if n_ids > 1
        ii = round(1 + (length(opts.cm) - 1) .* ...
            linspace(0.15, 0.85, n_ids));
        opts.cm = flip(opts.cm(ii,:), 1);
    else
        opts.cm = opts.cm(end - round(opts.cloc * size(opts.cm,1)), :);
    end
end

% Sphere geometry
[X, Y, Z] = sphere(60);

% Plot each primary particle
h_pp = gobjects(n_pp, 1);
for i = 1:n_pp
    h_pp(i) = surf(X * d_pp(i)/2 + x_pp(i), ...
                   Y * d_pp(i)/2 + y_pp(i), ...
                   Z * d_pp(i)/2 + z_pp(i));

    h_pp(i).EdgeColor = 'none';
    h_pp(i).FaceAlpha = opts.ft;
    h_pp(i).FaceLighting = 'gouraud';
    h_pp(i).AmbientStrength = 0.8;
    h_pp(i).DiffuseStrength = 0.2;
    h_pp(i).SpecularStrength = 0.05;
    h_pp(i).SpecularExponent = 2;
    h_pp(i).BackFaceLighting = 'lit';

    % Apply color
    if use_scalar_coloring
        h_pp(i).CData = ones(size(Z)) * colorval(i);  % map scalar across surface
        h_pp(i).CDataMapping = 'scaled';
        h_pp(i).FaceColor = 'flat';
    else
        h_pp(i).FaceColor = opts.cm(ids == colorval(i), :);
    end

    hold on
end

% General appearance
axis equal
axis off
grid off
camlight('right')

% Clear output if not needed
if nargout == 0
    clear h_pp
end
end
