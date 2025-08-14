function h_pp = PLOTPP_CONTINUOUS(x_pp, y_pp, z_pp, d_pp, colorval, opts)
% PLOTPP_CONTINUOUS displays primary particles in 3D with scalar or group coloring.
% ------------------------------------------------------------------------
% Inputs:
%   x_pp, y_pp, z_pp : Coordinates of primary particles
%   d_pp             : Diameters of primary particles
%   colorval         : Scalar (e.g., shielding) or group IDs for coloring
%   opts             : Struct with fields:
%                      - cm   : custom colormap (default = gray)
%                      - ft   : face transparency (default = 1)
%                      - cc   : 'on' for scalar colormap coloring
%                      - cloc : location for selecting default colormap color
%                      - clim : [min max] color limits (optional, for scalar)
% ------------------------------------------------------------------------
% Output:
%   h_pp             : Handles to the rendered particle surfaces
% ------------------------------------------------------------------------

% Number of primary particles
n_pp = numel(x_pp);

% --- Handle missing options ---
if ~exist('opts', 'var'), opts = struct(); end
if ~isfield(opts, 'cm'), opts.cm = gray; end
if ~isfield(opts, 'ft'), opts.ft = 1; end
if ~isfield(opts, 'cc'), opts.cc = 'off'; end
if ~isfield(opts, 'cloc'), opts.cloc = 0.25; end

% --- Determine coloring mode ---
% use_scalar_coloring = strcmpi(opts.cc, 'on') && isnumeric(colorval) && ...
%                       (any(colorval ~= round(colorval)) || any(colorval < 1));
use_scalar_coloring = strcmpi(opts.cc, 'on');

% --- Prepare color mapping ---
if use_scalar_coloring
    % Determine color limits (either user-specified or data-driven)
    if isfield(opts, 'clim') && numel(opts.clim) == 2
        vmin = opts.clim(1);
        vmax = opts.clim(2);
    else
        vmin = min(colorval);
        vmax = max(colorval);
    end
    cdata_mapping = 'scaled';
else
    % Discrete group-based coloring
    ids = unique(colorval);
    n_ids = length(ids);
    if n_ids > 1
        % Subsample and flip colormap range
        ii = round(1 + (length(opts.cm) - 1) .* linspace(0.15, 0.85, n_ids));
        opts.cm = flip(opts.cm(ii,:), 1);
    else
        % Single color fallback
        opts.cm = opts.cm(end - round(opts.cloc * size(opts.cm,1)), :);
    end
    cdata_mapping = 'direct';
end

% --- Generate sphere mesh for plotting ---
[X, Y, Z] = sphere(60);

% --- Initialize output handle array ---
h_pp = gobjects(n_pp, 1);

% --- Main rendering loop ---
for i = 1:n_pp
    % Generate surface of the i-th primary particle
    h_pp(i) = surf(X * d_pp(i)/2 + x_pp(i), ...
                   Y * d_pp(i)/2 + y_pp(i), ...
                   Z * d_pp(i)/2 + z_pp(i));

    % Surface appearance
    h_pp(i).EdgeColor = 'none';
    h_pp(i).FaceAlpha = opts.ft;
    h_pp(i).FaceLighting = 'gouraud';
    h_pp(i).AmbientStrength = 0.8;
    h_pp(i).DiffuseStrength = 0.2;
    h_pp(i).SpecularStrength = 0.05;
    h_pp(i).SpecularExponent = 2;
    h_pp(i).BackFaceLighting = 'lit';

    % --- Apply color mapping ---
    if use_scalar_coloring
        % Scalar-based coloring: assign value and scale
        h_pp(i).CData = ones(size(Z)) * colorval(i);
        h_pp(i).FaceColor = 'flat';
        h_pp(i).CDataMapping = 'scaled';

        % Make sure current axis reflects correct colormap
        ax = ancestor(h_pp(i), 'axes');
        colormap(ax, opts.cm)
        caxis(ax, [vmin vmax])
    else
        % Group-based coloring: direct color lookup
        h_pp(i).FaceColor = opts.cm(ids == colorval(i), :);
    end

    hold on
end

% --- Final appearance settings ---
axis equal
axis off
grid off
camlight('right')

% Clear unused handle if not requested
if nargout == 0
    clear h_pp
end
end