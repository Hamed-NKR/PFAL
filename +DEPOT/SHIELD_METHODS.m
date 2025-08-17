% the final calculation is wrong...

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
    
    % process in chunks to reduce memory usage
    ind_chunk = 1 : opts.chunkSize : length(kkk);
    dr_chunk = cell(numel(ind_chunk),1);    
    for c = ind_chunk

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
        dr_chunk{c} = r2 < r_thresh2;

        % count the number of screened circumferential points
        spp0(c:c_end) = sum(dr_chunk{c}, 2);
        
    end
    
    % concatinate dr_chunk
    dr_chunk = cat(1,dr_chunk{:});

    % count how many primary particles are along each discretization point
    n0_points = sum(dr_chunk, 1);

    % count screened points with the more rigorous layered criteria (i.e.
        % ...considers screening when at least two primary particles fall...
        % ...along discretization point)
    spp00 = sum(dr_chunk(:, n0_points>1), 2);

    % select maximum obscured perimeter points as final screening
    if ~isempty(spp0(dz(kkk)))
        spp_i_singleSide(k) = max(spp0(dz(kkk))) / n_points; % choose...
            % ...only the ones with lower z
    end

    spp_i_doubleSide(k) = max(spp0) / n_points; % choose regardless of z

    spp_i_doubleLayer(k) = max(spp00) / n_points; % choose considering...
        % ...a simple optical depth criterion

end

end

