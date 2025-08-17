function spp = SHIELD(pp, n_sample, method, opts)
% "SHIELD" computes screening factor of all individual primary particles...
%   ...in one aggregate under a specified screening rule from three...
%   ...available choices.
%
% ----------------------------------------------------------------------- %
% Inputs:
%   Required:
%       pp: [N x 6] with columns [id, d, x, y, z, subagg]
%   Optional (defaults if omitted/empty):
%       n_sample: perimeter samples per primary particle (default 64)
%       method: 'opaque' | 'transparent' | 'semitransparent'...
%           [default: 'semitransparent']
%       opts: struct for optimizing memory usage
%           [fields with defaults: .chunkSize (512), .xyPrune (true)]
%
% Screening rules:
%       'opaque': front-only (+z), overlay threshold >=1
%       'transparent': ignores z, overlay threshold >=1
%       'semitransparent': ignores z, overlay threshold >=2
%
% ----------------------------------------------------------------------- %
% Outputs:
%   spp: [N x 1] vector of screening factor for individual primary...
%       ...particles (0: visible, 1:fully covered)
%
% ----------------------------------------------------------------------- %

%%% setting defaults

% number of points on each primary particle's xy perimeter projection
if nargin < 2 || isempty(n_sample), n_sample = 64; end

% screening rule
if nargin < 3 || isempty(method), method = 'semitransparent'; end

% memory usage
def.chunkSize = 512; % break matrix calculations into chunks
def.xyPrune = true; % ignore calculation if primary particles...
    % ...do not overlay
if nargin < 4 || isempty(opts), opts = def; else
    if ~isfield(opts,'chunkSize') || isempty(opts.chunkSize)
        opts.chunkSize = def.chunkSize;
    end
    if ~isfield(opts,'xyPrune') || isempty(opts.xyPrune)
        opts.xyPrune   = def.xyPrune;
    end
end

npp = size(pp,1); % number of primary particles
spp = zeros(npp,1,'double'); % screening factor for each primary particle

d = pp(:,2); % primary particle diameter
r = 0.5*d; % primary particle radius

% primary particle center coordinates [x y z]
xc = pp(:,3); yc = pp(:,4); zc = pp(:,5);

% angles for perimeter sampling (uniform spacing)
ang = (0:n_sample-1) * (2*pi / max(1,n_sample));

% tolerances (scale with geometry)
tolr = 1e3*eps(max(r)); % ~1e-21 m for r~1e-8 m
tolz = 1e3*eps(max(abs(zc)) + max(r));

% tolz = 1e-14 * max(1, max(abs(zc)) + max(r)); % z direction
% tolr = 1e-14 * max(1, max(r)); % radial direction

for k = 1:npp % iteration over individual primary particles

    % generate perimeter points
    rk = r(k);
    xk = xc(k) + rk * cos(ang);
    yk = yc(k) + rk * sin(ang);
    
    % indices of "other" primary particles
    if npp <= 1, spp(k) = 0; continue; end
    tgt_all = setdiff(1:npp, k);

    % keep only "other" primary particles that overlay in xy plane with...
        % ...selected primary particle 
    if opts.xyPrune && ~isempty(tgt_all)
        dx0  = xc(tgt_all) - xc(k);
        dy0  = yc(tgt_all) - yc(k);
        rtgt = r(tgt_all);
        keep = (dx0.^2 + dy0.^2) <= (rk + rtgt + tolr).^2;
        tgt_all = tgt_all(keep);
        if isempty(tgt_all), spp(k) = 0; continue; end
    end

    % initialize per-point overlay counters
    switch method
        case 'opaque' % only +z occluders
            n_plus = zeros(1, n_sample, 'uint16'); 
        case 'transparent' % both +z and -z occluders
            n_all  = zeros(1, n_sample, 'uint16'); 
        case 'semitransparent' % minimum two occluders along z
            n_all  = zeros(1, n_sample, 'uint16');
        otherwise
            error('Unknown method: %s', method);
    end

    % stream pair-wise primary particle distance calculations in chunks...
        % ...for memory safety
    cs = min(opts.chunkSize, max(1, numel(tgt_all)));
    for c = 1:cs:numel(tgt_all)
        idx = tgt_all(c : min(c+cs-1, numel(tgt_all)));

        xt = xc(idx);
        yt = yc(idx);
        zt = zc(idx);
        rt = r(idx);

        % whether each perimeter points is inside "other" primary particles
        dx = xk - xt; % [numel(idx) x n_sample]
        dy = yk - yt; % [numel(idx) x n_sample]
        inside = (dx.^2 + dy.^2) < ((rt + tolr).^2);
        
        % count how many "other" primary particles overlay with each...
            % ...perimeter point
        switch method
            case 'opaque'
                front = (zt - zc(k)) > tolz;
                if any(front)
                    n_plus = n_plus + uint16(sum(inside(front,:), 1));
                end
            case 'transparent'
                n_all = n_all + uint16(sum(inside, 1));
            case 'semitransparent'
                n_all = n_all + uint16(sum(inside, 1));
        end

    end
    
    % calculate fraction of perimeter points that are occluded
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
