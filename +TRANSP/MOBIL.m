function pars = MOBIL(pars, fl, params_const, opts)
% "MOBIL" calculates the transport properties of aggregates.
% ----------------------------------------------------------------------- %
% 
% Inputs/output:
%     pars: Aggregates information structure containting their...
%         ...geometrical and diffusive properties
%     fl: Fluid info structure
%     params_const: Problem's table of constant physical properties
%     opts: Mobility calculation options 
% ----------------------------------------------------------------------- %

% initialize calculation method variable
if ~exist('opts', 'var') 
    opts = struct();
end
if (~isfield(opts, 'mtd')) || isempty(opts.mtd)
    opts.mtd = 'continuum'; % use continuum regime approximation
end

% initialize time-step sensitivity variable
if (~isfield(opts, 'c_dt')) || isempty(opts.c_dt)
    opts.c_dt = 1; % multiply random-walk timestep by this factor
end

% Total number of (independent) particles
if isa(pars, 'AGG')
    n_par = length(pars);
else
    n_par = size(pars.n, 1);
end

% Calculating mobility properties

% Mobility size
if strcmp(opts.mtd, 'continuum')
    npp = cat(1, pars.n);
    dm = 1.29 * npp.^(-0.13) .* cat(1, pars.dg);
elseif strcmp(opts.mtd, 'interp')
    opts_proj.tbar = 'off';
    da = 2 * sqrt(PAR.PROJECTION(pars, [], 1e2, 5, [], opts_proj) / pi);
    dm = TRANSP.DIAMOBIL(pars.dg, da, fl);
end

% Cunningham correction factor (-)
kn_kin = (2 * fl.lambda) ./ dm; % Kinetic (momentum) Knudsen number (-)
alpha = 1.254; 
beta = 0.4;
gamma = 1.1;
cc = 1 + kn_kin .* (alpha + beta .* exp(-gamma ./ kn_kin));
    

% Aggregate mass (kg)
m = ones(n_par,1);
for i = 1 : n_par
    if isa(pars, 'AGG')
        m(i) = (params_const.Value(1) * pi / 6) * sum((pars(i).pp.d(:)).^3);
    else
        m(i) = (params_const.Value(1) * pi / 6) * sum((pars.pp{i}(:,2)).^3);
    end
end

rho = (6 / pi) * (m ./ (dm.^3)); % Effective density (kg/m3)
tau = rho .* (dm.^2) .* cc ./ (18 * fl.mu); % Response (relaxation) time (s)
kb = params_const.Value(3); % Boltzmann's constant (j/k)
f = (3 * pi * fl.mu) .* dm ./ cc; % Friction factor (-)
delt = opts.c_dt * f .* (dm.^2) ./ (6 * kb * fl.temp); % Marching timestep
diff = (kb * (fl.temp)) ./ f; % Particle diffusivity (m2/s)

lambda = sqrt(m .* kb .* (fl.temp)) ./ f; % Diffusive mean free path (m)
kn_diff = 2 * lambda ./ dm; % Diffusive Knudsen number (-)

% Updating the transport properties for the structure/class of aggregates
if isa(pars, 'AGG')
    for i = 1 : n_par
        pars(i).rho = rho(i);
        pars(i).m = m(i);
        pars(i).tau = tau(i);
        pars(i).f = f(i);
        pars(i).delt = delt(i);
        pars(i).diff = diff(i);
        pars(i).lambda = lambda(i);
        pars(i).kn_kin = kn_kin(i);
        pars(i).kn_diff = kn_diff(i);
    end

else
    pars.rho = rho;
    pars.m = m;
    pars.dm = dm;
    if strcmp(opts.mtd, 'interp'); pars.da = da; end
    pars.tau = tau;
    pars.f = f;
    pars.delt = delt;
    pars.diff = diff;
    pars.lambda = lambda;
    pars.kn_kin = kn_kin;
    pars.kn_diff = kn_diff;
end

end
