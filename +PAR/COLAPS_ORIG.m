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
defSpecs = struct('k_spring', 0.2, 'k_decay', 1.0, 'jit', 0.0,...
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

