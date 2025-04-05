function h_pp = PLOTPP_CONTINUOUS(x_pp, y_pp, z_pp, d_pp, colorval, opts)
% "PLOTPP" displays the 3d schematic of primary particles in 3d space.
% ----------------------------------------------------------------------- %
% Input:
%     x_pp: x coordinate of primaries 
%     y_pp: y ~
%     z_pp: z ~
%     d_pp: Diameter of primaries
%     colorval: either group ids or scalar value (e.g. shielding)
%     opts: plotting options
% ----------------------------------------------------------------------- %
% Output:
%     h_pp: Figure handle
% ----------------------------------------------------------------------- %

n_pp = numel(x_pp); % Number of primary particles

if ~exist('opts', 'var'); opts = struct(); end
if ~isfield(opts, 'cm'); opts.cm = gray; end
if ~isfield(opts, 'ft'); opts.ft = 1; end
if ~isfield(opts, 'cc'); opts.cc = 'off'; end
if ~isfield(opts, 'cloc'); opts.cloc = 0.25; end

% Use scalar mode if input has continuous values
use_scalar_coloring = strcmpi(opts.cc, 'on') && isnumeric(colorval) && ...
                      (any(colorval ~= round(colorval)) || any(colorval < 1));

if use_scalar_coloring
    % Normalize color values to [0, 1]
    vmin = min(colorval);
    vmax = max(colorval);
    if vmax == vmin
        colorval_norm = ones(size(colorval));  % avoid divide-by-zero
    else
        colorval_norm = (colorval - vmin) / (vmax - vmin);
    end
    cm = opts.cm;
else
    ids = unique(colorval);
    n_agg = length(ids);
    if n_agg > 1
        ii = round(1 + (length(opts.cm) - 1) .* ...
                  (0.15 : 0.7 / (n_agg - 1) : 0.85)');
        cm = flip(opts.cm(ii,:),1);
    else
        cm = opts.cm(end - round(opts.cloc * length(opts.cm)),:);
    end
end

[X, Y, Z] = sphere(60);

for i = 1 : n_pp
    h_pp = surf(X .* d_pp(i)/2 + x_pp(i), ...
                Y .* d_pp(i)/2 + y_pp(i), ...
                Z .* d_pp(i)/2 + z_pp(i)); % plot primaries

    h_pp.EdgeColor = 'none';
    h_pp.FaceAlpha = opts.ft;
    h_pp.FaceLighting = 'gouraud';
    h_pp.AmbientStrength = 0.8;
    h_pp.DiffuseStrength = 0.2;
    h_pp.SpecularStrength = 0.05;
    h_pp.SpecularExponent = 2;
    h_pp.BackFaceLighting = 'lit';

    % Assign color
    if use_scalar_coloring
        cidx = max(1, round(colorval_norm(i) * (size(cm,1)-1))) + 1;
        h_pp.FaceColor = cm(cidx,:);
    else
        h_pp.FaceColor = cm(ids == colorval(i),:);
    end

    hold on
end

axis equal
grid off
axis off
camlight('right');

if nargout == 0
    clear h_pp;
end
end
