function h_pp = PLOTPP(x_pp, y_pp, z_pp, d_pp, id_agg, opts)
% "PLOTPP" displays the 3d schematic of primary particles in 3d space.
% ----------------------------------------------------------------------- %
% 
% Input:
%     x_pp: x coordinate of primaries 
%     y_pp: y ~
%     z_pp: z ~
%     d_pp: Dimater of primaries
%     id_agg: ids of primaries (which aggregate they belong to)
%     opts: plotting options
% ----------------------------------------------------------------------- %
% 
% Output:
%     h_pp: The output figure handle
% ----------------------------------------------------------------------- %

n_pp = numel(x_pp); % Number of primary particles

% Set variable for plot properties
if ~exist('opts', 'var') 
    opts = struct();
end

% Setting default particle colormap
if ~isfield(opts, 'cm'); opts.cm = []; end
if isempty(opts.cm); opts.cm = gray; end

% Setting default particle transparency as 1 (fully opaque)
if ~isfield(opts, 'ft'); opts.ft = []; end
if isempty(opts.ft); opts.ft = 1; end

% Setting defaults on whether colorcode the aggs
if ~isfield(opts, 'cc'); opts.cc = 'off'; end
if ismember(opts.cc, {'OFF', 'Off', 'off'}) || ~exist('id_agg', 'var') ||...
        isempty(id_agg)
    id_agg = ones(n_pp,1);
end

% Setting the location of colormaps to choose for painting a uniform...
% ...aggregates
if ~isfield(opts, 'cloc'); opts.cloc = []; end
if isempty(opts.cloc); opts.cloc = 0.1; end

ids = unique(id_agg); % unique agg ids existing in the population
n_agg = length(ids); % number of aggs existing

% make colormap
if n_agg > 1
    ii = round(1 + (length(opts.cm) - 1) .* (0.15 : 0.7 / (n_agg - 1) : 0.85)');
    cm = opts.cm(ii,:);
    cm = flip(cm,1);
else
    cm = opts.cm(end - round(opts.cloc * length(opts.cm)),:);
end

% Plotting aggregates
[X,Y,Z] = sphere(60);

for i = 1 : n_pp
    h_pp = surf(X .* d_pp(i) ./ 2 + x_pp(i),...
        Y .* d_pp(i) ./ 2 + y_pp(i), Z .* d_pp(i) ./ 2 + z_pp(i)); % plot primaries
    
    % set graphics
    h_pp.EdgeColor = 'none';
    if n_agg > 1
        h_pp.FaceColor = cm(ids == id_agg(i),:);
    else
        h_pp.FaceColor = cm;
    end
    h_pp.FaceAlpha = opts.ft;
    h_pp.FaceLighting = 'gouraud';
    h_pp.AmbientStrength = 0.8;
    h_pp.DiffuseStrength = 0.2;
    h_pp.SpecularStrength = 0.05;
    h_pp.SpecularExponent = 2;
    h_pp.BackFaceLighting = 'lit';
    
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

