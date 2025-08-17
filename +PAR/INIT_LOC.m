function [pars, params_ud] = INIT_LOC(pars, params_ud, n_itr, opts)
% "INIT_LOC" randomly distributes the particles throughout the...
%     ...computational domain.
% ----------------------------------------------------------------------- %
% 
% Inputs/Output:
%     pars: Particle information structure
%     params_ud: User defined parameters (including the domain size)
%     n_itr: maximum number of iterations to locate the particles
% ----------------------------------------------------------------------- %

% Compiling pp info
if isa(pars, 'AGG')
    pp = AGG.COMPILEPP(pars);
else
    pp = cell2mat(pars.pp);
end

% Initializing number of iterations if not defined
if ~exist('n_itr', 'var') || isempty(n_itr)
    n_itr = 1e4;
end

n_par = size(pars.n,1); % Total number of (independent) particles
pars.r = PAR.COM(pars.pp, pars.n); % Assigning center of mass as initial...
    % ...location of the aggregates
dmax = PAR.TERRITORY(pars.pp, pars.n); % Maximum distance from the...
    % ...center of each aggregate

% option for adjusting volume fraction
if ~exist('opts', 'var') 
    opts = struct();
end
if (~isfield(opts, 'vf')) || isempty(opts.vf)
    opts.vf = 'off'; % disable domain size adjusting
end
if (~isfield(opts, 'vf_type')) || isempty(opts.vf_type)
    opts.vf_type = 'pp'; % set the domain size adjusting to be based on...
        % ...volume fraction of primary particles
end

% Updating the domain size based on the volume fraction 
if (strcmp(opts.vf, 'on') || strcmp(opts.vf, 'On') ||...
        strcmp(opts.vf, 'ON')) && isa(params_ud.Value(1), 'double') &&...
        (params_ud.Value(1) > 0)
    if strcmp(opts.vf_type, 'pp')
        params_ud.Value(2) = (pi * sum(pp(:,2).^3) /...
            (6 * params_ud.Value(1)))^(1/3);
    else
        params_ud.Value(2) = (pi * sum(dmax.^3) /...
            (6 * params_ud.Value(1)))^(1/3);
    end
    params_ud.Value(3) = params_ud.Value(2);
    params_ud.Value(4) = params_ud.Value(2);
end

% Initialization of the location array
pars.r = rand(n_par,3) .* (repmat((params_ud.Value(2:4))',n_par,1) -...
     repmat(dmax,1,3))+ (repmat(dmax,1,3) ./ 2);

% Making particle pair indices
ind_pars = (1:n_par)';
ind_pars = [repelem(ind_pars,n_par,1), repmat(ind_pars,n_par,1)];

% identifying repeating pairs
rmv1 = (1:n_par)';
rmv1 = repelem((rmv1-1).*n_par,1:n_par);
rmv2 = repmat((1:n_par)',[1 n_par]);
rmv2 = triu(rmv2);
rmv2 = reshape(rmv2,n_par^2,1);
rmv2(rmv2 == 0) = [];
rmv = rmv1 + rmv2;
ind_pars(rmv,:) = [];

% Generating the "OVR" inputs:
d_pars = [repelem(dmax, n_par, 1), repmat(dmax, n_par, 1)];
    % size input
d_pars(rmv,:) = [];
r_pars = [repelem(pars.r, n_par, 1), repmat(pars.r, n_par, 1)];
    % location input
r_pars(rmv,:) = [];

ovrs = COL.OVR(r_pars, d_pars); % Checking initial overlapping...
% ...between the particles
ind_err = 0; % Initializing error generation index

fprintf('Initializing particle locations based on maximum extents...')
disp(' ')
UTILS.TEXTBAR([0, n_itr]); % Initialize textbar
UTILS.TEXTBAR([1, n_itr]); % Iteration 1 already done

% Reinitializing overlapped particles
while ~ isempty(find(ovrs == 1, 1))
    if ind_err >= n_itr
        error('Failed locating particles...')
    end
    
    % Updating the location of overlapped particles
    ind_updt = ind_pars(ovrs == 1, :); % Indices of updated particles
    ind_updt = unique(ind_updt(:)); % removing repeating indices
    pars.r(ind_updt, 1:3) = rand(size(ind_updt,1),3) .*...
        (repmat((params_ud.Value(2:4))', size(ind_updt,1),1) -...
        repmat(dmax(ind_updt),1,3)) +...
        (repmat(dmax(ind_updt),1,3) ./ 2);
    r_pars = [repelem(pars.r,n_par,1), repmat(pars.r,n_par,1)];
    r_pars(rmv,:) = [];
    ovrs = COL.OVR(r_pars, d_pars); % Rechecking the overlapping
    
    ind_err = ind_err + 1; % Updating error index
    
    UTILS.TEXTBAR([ind_err, n_itr]); % Update textbar
end

disp(' ')

% Updating the primary particle locations based on their new random...
    % ...center positions
pp_rc = PAR.COM(pars.pp, pars.n); % Initial particle centers
dr0 = pars.r - pp_rc; % Translation vectors
pars.pp = PAR.TRANSLATE(pars.pp, pars.r, pars.n, dr0);

end
