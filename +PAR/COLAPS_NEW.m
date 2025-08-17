function [pps, n_steps] = COLAPS(pp0, n_steps, specs, opts)
% Minimal, fast collapse: global spring (decays over time) + LJ with smooth cutoff.
% Accepts small overlaps; tuned to compact in ~1e4–2e4 steps.

% ------------------ Defaults ------------------
defSpecs = struct( ...
    'k_spring0', 0.25, ...   % initial spring
    'k_springf', 0.08,  ...  % final spring at last step
    'jit', 0.0, 'jit_decay', 0.98, ...
    'lj_eps',   0.80, ...    % cohesion strength (lower = easier compaction)
    'lj_decay', 0.94, ...    % damping (friction)
    'dt',       0.004, ...   % time step
    'rc_fac',   1.6,  ...    % smooth-taper start: rc = rc_fac * sigma
    'rcut_fac', 2.2   ...    % smooth-taper end:  rcut = rcut_fac * sigma
);
defOpts  = struct('showProgress', true);

if nargin < 2 || isempty(n_steps), n_steps = 2e4; end
if nargin < 3 || isempty(specs),   specs = defSpecs; else
    fn = fieldnames(defSpecs); for i=1:numel(fn)
        if ~isfield(specs,fn{i}), specs.(fn{i}) = defSpecs.(fn{i}); end
    end
end
if nargin < 4 || isempty(opts),    opts  = defOpts; else
    fn = fieldnames(defOpts); for i=1:numel(fn)
        if ~isfield(opts,fn{i}),  opts.(fn{i})  = defOpts.(fn{i}); end
    end
end

% ------------------ Normalize inputs ------------------
d0  = pp0(:,2);             % diameters
x0  = pp0(:,3:5);           % coordinates
xcm0 = sum((d0.^3).*x0, 1) ./ sum(d0.^3);
r0   = d0/2; rm0 = mean(r0);

r = r0 / rm0;               % normalized radii
x = x0 / rm0; x = x - mean(x,1);

npp = size(x,1);
v   = zeros(npp,3);
m   = (4/3)*pi*(r.^3);

% pair params that never change
sig = (r + r.') / (2^(1/6));                 % LJ minimum at contact
rc   = specs.rc_fac  .* sig;
rcut = specs.rcut_fac.* sig;

pps = cell(n_steps,1);

if opts.showProgress
    disp('Aggregate collapsing (simple mode)...');
    UTILS.TEXTBAR([0, n_steps]);
end

% ------------------ Time stepping ------------------
for kk = 1:n_steps

    % Exponential schedule for spring strength (simple, robust)
    if specs.k_spring0 > 0
        t = kk / n_steps;
        k_spring = specs.k_spring0 * (specs.k_springf/specs.k_spring0)^t;
    else
        k_spring = 0;
    end

    % jitter decay (usually 0 while debugging)
    specs.jit = specs.jit * specs.jit_decay;

    % Pairwise geometry
    dvec = permute(x,[1 3 2]) - permute(x,[3 1 2]);   % (i,j,xyz)
    ds   = sqrt(sum(dvec.^2,3));                      % distances
    ds(1:npp+1:end) = Inf;
    d_unit = dvec ./ ds;

    % Forces: spring-to-COM (mass-weighted, polydisperse)
    cen = sum(x .* (r.^3), 1) ./ sum(r.^3);
    F   = -k_spring * (x - cen);

    % Lennard–Jones with smooth taper (single pass)
    eps0   = 1e-12;
    Fvwmag = 48*specs.lj_eps .* ( (sig.^12) ./ (ds.^13 + eps0) ...
                                - 0.5*(sig.^6) ./ (ds.^7  + eps0) );

    % C2 smooth taper rc→rcut (avoids force jumps & helps reconnection)
    S = ones(size(ds));
    mask = ds > rc & ds < rcut;
    xi = zeros(size(ds));
    xi(mask) = (ds(mask) - rc(mask)) ./ (rcut(mask) - rc(mask)); % 0..1
    S(mask)  = 1 - 3*xi(mask).^2 + 2*xi(mask).^3;
    S(ds >= rcut) = 0;
    Fvwmag = Fvwmag .* S;

    % add vdW once
    F = F + squeeze(sum(repmat(Fvwmag,1,1,3) .* d_unit, 2)) ...
          + specs.jit * randn(npp,3);

    % Leapfrog
    a = F ./ m;
    v = v + a * specs.dt;
    x = x + v * specs.dt;
    v = v * specs.lj_decay;

    % Un-normalize for output
    xt   = rm0 * x;
    xcmT = sum((d0.^3).*xt, 1) ./ sum(d0.^3);
    xt   = xt + (xcm0 - xcmT);

    pps{kk} = pp0;
    pps{kk}(:,3:5) = xt;

    if opts.showProgress, UTILS.TEXTBAR([kk, n_steps]); end
end
end
