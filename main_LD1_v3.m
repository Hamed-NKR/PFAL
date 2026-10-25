clc
clear
clf('reset')
close all

%% Load first-stage LD configuration %%

% First-stage LD now reads all user-tuned parameters from JSON so the
% aggregate-library run can be reproduced without editing this script.
cfg_ld1 = UTILS.LOAD_MAIN_LD1_CONFIG;
cfg_results = cfg_ld1.results;
cfg_sim = cfg_ld1.simulation;
cfg_projection = cfg_ld1.projection;
cfg_transport = cfg_ld1.transport;

if ~isempty(cfg_sim.rng_seed)
    rng(cfg_sim.rng_seed);
end

[params0_ud, params0_const] = UTILS.LD1_PARAMS_FROM_CONFIG( ...
    cfg_transport.domain, cfg_transport.particles, cfg_transport.fluid);

% Initialize particle and fluid structures using the same table contract as
% the legacy PFA first-stage workflow.
[pars0, fl0] = TRANSP.INIT_DOM(params0_ud, params0_const);

opts0_fl = cfg_transport.fluid_options;
[fl0.mu, fl0.lambda] = TRANSP.FLPROPS(fl0, params0_const, opts0_fl);

%% Initialize simulation variables %%

n_temporal = cfg_sim.n_temporal; % number of saved growth states per trial
n_trial = cfg_sim.n_trial; % number of independent LD trials
s_pp0 = cfg_sim.primary_size_gsd_between_trials; % between-trial pp GM spread
npp0_min = cfg_sim.npp_min; % smallest aggregate size target to store
npp0_max = cfg_sim.npp_max; % largest aggregate size target to store
j_max = cfg_sim.j_max; % maximum number of iterations per trial

opts0_mobil = cfg_transport.mobility_options;
opts0_grow = cfg_transport.growth_options;
opts0_loc = cfg_transport.location_options;
opts0_prj = cfg_projection.options;

n0_mc_prj = cfg_projection.n_mc;
n0_ang_prj = cfg_projection.n_ang;

% Make a 2D cell array to store primary-particle data over time and trials.
pp0 = cell(n_temporal, n_trial);
pp0_n = cell(n_temporal, n_trial);

% Allocate aggregate-property snapshots aligned with pp0.
parsdata0 = struct('dpp', cell(n_temporal, n_trial), ...
    'sigmapp', cell(n_temporal, n_trial), 'da', cell(n_temporal, n_trial), ...
    'dm', cell(n_temporal, n_trial), 'dg', cell(n_temporal, n_trial));

% Aggregate-size targets for intermediate storage.
r_temporal = (npp0_max / npp0_min)^(1 / (n_temporal - 1));
npp0_temporal = zeros(n_temporal, 1);
for i = 1 : n_temporal
    npp0_temporal(i) = round(npp0_min * (r_temporal^(i - 1)));
end

fprintf('Initializing first-stage LD aggregate-library run...\n')
fprintf('LD1 config: %s\n', cfg_ld1.config_file)

% Generate a lognormal distribution for the ensemble geometric mean of
% primary-particle diameter in each independent trial.
dpp0_ens = lognrnd(log(params0_ud.Value(8)), log(s_pp0), [n_trial, 1]);

% Allocate ensemble-averaged properties at every timestep.
ensdata0.n_agg = zeros(j_max, n_trial);
ensdata0.t = zeros(j_max, n_trial);
ensdata0.tau = zeros(j_max, 2 * n_trial);
ensdata0.kn_kin = zeros(j_max, 2 * n_trial);
ensdata0.kn_diff = zeros(j_max, 2 * n_trial);
ensdata0.dpp = zeros(j_max, 2 * n_trial);
ensdata0.sigmapp = zeros(j_max, 2 * n_trial);
ensdata0.da = zeros(j_max, 2 * n_trial);
ensdata0.dm = zeros(j_max, 2 * n_trial);

%% Perform first-stage Langevin dynamics simulations %%

for i = 1 : n_trial
    fprintf('trial %d:\n', i)

    % Initialize primary-particle diameters. This mirrors the first-stage
    % setup in main_PFA_v2_2.m while allowing per-run values from JSON.
    [pp0_d, pars0.n] = PAR.INIT_DIAM(params0_ud.Value(5), ...
        params0_ud.Value(6:7), [dpp0_ens(i); params0_ud.Value(9:10)]);

    pars0.pp = mat2cell([(1:size(pp0_d, 1))', pp0_d, ...
        zeros(size(pp0_d, 1), 3), (1:size(pp0_d, 1))'], pars0.n);

    if params0_ud.Value(6) ~= 0
        pars0.pp = PAR.INIT_MORPH_RAND(pars0.pp);
    end

    [pars0, params0_ud] = PAR.INIT_LOC(pars0, params0_ud, [], opts0_loc);
    fl0.size = params0_ud.Value(2:4);

    pars0 = PAR.SIZING(pars0);
    pars0 = TRANSP.MOBIL(pars0, fl0, params0_const, opts0_mobil);
    ensdata0 = record_ensemble(ensdata0, 1, i, pars0, opts0_mobil);

    pars0.v = PAR.INIT_VEL(pars0.pp, pars0.n, fl0.temp, params0_const);

    j = 2;
    jj = 1;
    seen_aggregates = containers.Map('KeyType', 'char', 'ValueType', 'logical');

    fprintf('Simulating...\n')
    UTILS.TEXTBAR([0, j_max]);
    UTILS.TEXTBAR([1, j_max]);

    while (j <= j_max) && ...
            (sum(cat(1, pars0.n) >= npp0_max) < ...
            round(0.9 * length(cat(1, pars0.n)))) && ...
            (length(pars0.n) > 1)

        pars0 = TRANSP.MARCH(pars0, fl0, params0_const);
        pars0 = TRANSP.PBC(params0_ud.Value(2:4), pars0);
        pars0 = COL.GROW(pars0, opts0_grow);
        pars0 = PAR.SIZING(pars0);
        pars0 = TRANSP.MOBIL(pars0, fl0, params0_const, opts0_mobil);

        if (jj <= n_temporal) && ...
                (sum(cat(1, pars0.n) >= npp0_temporal(jj)) >= ...
                round(0.9 * length(cat(1, pars0.n))))
            [pp0, pp0_n, parsdata0, seen_aggregates] = store_snapshot( ...
                pp0, pp0_n, parsdata0, seen_aggregates, pars0, ...
                jj, i, n_temporal, params0_ud.Value(5), opts0_mobil, ...
                n0_mc_prj, n0_ang_prj, opts0_prj);
            jj = jj + 1;
        end

        ensdata0 = record_ensemble(ensdata0, j, i, pars0, opts0_mobil);

        UTILS.TEXTBAR([j, j_max]);
        j = j + 1;
    end

    if jj <= n_temporal
        [pp0, pp0_n, parsdata0, seen_aggregates] = store_snapshot( ...
            pp0, pp0_n, parsdata0, seen_aggregates, pars0, ...
            jj, i, n_temporal, params0_ud.Value(5), opts0_mobil, ...
            n0_mc_prj, n0_ang_prj, opts0_prj);
    end

    % Remove unused real-time rows for this trial only after the run is done.
    ensdata0 = trim_trial_ensemble(ensdata0, j, i);
    disp(newline)
end

%% Compile and save the aggregate library %%

pp0_nested = pp0;
pp0_n_nested = pp0_n;
parsdata0_nested = parsdata0;

pp0 = flatten_saved_pp(pp0_nested);
pp0_n = flatten_saved_counts(pp0_n_nested);

if isempty(pp0)
    error('PFAL:main_LD1_v3:EmptyLibrary', ...
        'LD1 produced no aggregate snapshots. Increase j_max or inspect the config.')
end

pars_ld1.pp = pp0;
pars_ld1.n = pp0_n;
pars_ld1 = PAR.SIZING(pars_ld1);

pars_ld1.da = flatten_parsdata_field(parsdata0_nested, 'da');
pars_ld1.dm = flatten_parsdata_field(parsdata0_nested, 'dm');
ensdata0 = trim_empty_ensemble_rows(ensdata0);

metadata = struct();
metadata.source_script = mfilename('fullpath');
metadata.config_file = cfg_ld1.config_file;
metadata.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
metadata.run_label = cfg_results.run_label;
metadata.output_file = fullfile(cfg_results.root, cfg_results.library_file);
metadata.n_aggregate = numel(pp0);
metadata.n_primary = sum(pp0_n);
metadata.primary_size_gsd_between_trials = s_pp0;
metadata.primary_particle_gsd_between_aggregates = ...
    cfg_transport.particles.d_pp_gsd_between_aggregates;
metadata.storage_targets_npp = npp0_temporal;

if ~isfolder(cfg_results.root)
    mkdir(cfg_results.root);
end

save(metadata.output_file, 'pp0', 'pp0_n', 'pars_ld1', 'ensdata0', ...
    'parsdata0', 'pp0_nested', 'pp0_n_nested', 'parsdata0_nested', ...
    'cfg_ld1', 'metadata', '-v7.3')

fprintf('Saved LD1 aggregate library: %s\n', metadata.output_file)

%% Local functions %%

function ensdata0 = record_ensemble(ensdata0, row_idx, trial_idx, pars0, opts0_mobil)
col = 2 * (trial_idx - 1) + (1:2);
ensdata0.t(row_idx, trial_idx) = min(pars0.delt);
if row_idx > 1
    ensdata0.t(row_idx, trial_idx) = ...
        ensdata0.t(row_idx - 1, trial_idx) + min(pars0.delt);
end
ensdata0.n_agg(row_idx, trial_idx) = length(pars0.n);
ensdata0.tau(row_idx, col) = [mean(pars0.tau), std(pars0.tau)];
ensdata0.kn_kin(row_idx, col) = [mean(pars0.kn_kin), std(pars0.kn_kin)];
ensdata0.kn_diff(row_idx, col) = [mean(pars0.kn_diff), std(pars0.kn_diff)];
ensdata0.dpp(row_idx, col) = [geomean(pars0.dpp_g(:, 1)), ...
    UTILS.GEOSTD(pars0.dpp_g(:, 1))];
ensdata0.sigmapp(row_idx, col) = [geomean(pars0.dpp_g(:, 2)), ...
    UTILS.GEOSTD(pars0.dpp_g(:, 2))];
ensdata0.dm(row_idx, col) = [geomean(pars0.dm), UTILS.GEOSTD(pars0.dm)];
if strcmp(opts0_mobil.mtd, 'interp') && isfield(pars0, 'da') && ~isempty(pars0.da)
    ensdata0.da(row_idx, col) = [geomean(pars0.da), UTILS.GEOSTD(pars0.da)];
end
end

function [pp0, pp0_n, parsdata0, seen_aggregates] = store_snapshot( ...
    pp0, pp0_n, parsdata0, seen_aggregates, pars0, storage_idx, trial_idx, ...
    n_temporal, n_initial, opts0_mobil, n0_mc_prj, n0_ang_prj, opts0_prj)

if ~strcmp(opts0_mobil.mtd, 'interp')
    pars0.da = 2 * sqrt(PAR.PROJECTION(pars0, [], n0_mc_prj, ...
        n0_ang_prj, [], opts0_prj) / pi);
end

keep = false(length(pars0.n), 1);
for k = 1 : length(pars0.n)
    signature = aggregate_signature(pars0.pp{k}, n_initial);
    if ~isKey(seen_aggregates, signature)
        seen_aggregates(signature) = true;
        keep(k) = true;
    end
end

if ~any(keep)
    return
end

pp_save = pars0.pp(keep);
n_save = pars0.n(keep);
id_offset = ((trial_idx - 1) * n_temporal + storage_idx - 1) * n_initial;
if id_offset ~= 0
    for k = 1 : numel(pp_save)
        pp_save{k}(:, 1) = pp_save{k}(:, 1) + id_offset;
    end
end

pp0{storage_idx, trial_idx} = pp_save;
pp0_n{storage_idx, trial_idx} = n_save;
parsdata0(storage_idx, trial_idx).dpp = pars0.dpp_g(keep, 1);
parsdata0(storage_idx, trial_idx).sigmapp = pars0.dpp_g(keep, 2);
parsdata0(storage_idx, trial_idx).da = pars0.da(keep);
parsdata0(storage_idx, trial_idx).dm = pars0.dm(keep);
parsdata0(storage_idx, trial_idx).dg = pars0.dg(keep);
end

function signature = aggregate_signature(pp, n_initial)
ids = mod(pp(:, 1), n_initial);
ids(ids == 0) = n_initial;
ids = sort(ids(:));
signature = sprintf('%d_', ids);
end

function pp_flat = flatten_saved_pp(pp_cells)
pp_cells = pp_cells(:);
pp_cells = pp_cells(~cellfun('isempty', pp_cells));
if isempty(pp_cells)
    pp_flat = {};
    return
end
pp_flat = cat(1, pp_cells{:});
pp_flat = pp_flat(:);
pp_flat = pp_flat(~cellfun('isempty', pp_flat));
end

function counts_flat = flatten_saved_counts(count_cells)
count_cells = count_cells(:);
count_cells = count_cells(~cellfun('isempty', count_cells));
if isempty(count_cells)
    counts_flat = zeros(0, 1);
else
    counts_flat = cat(1, count_cells{:});
end
end

function values = flatten_parsdata_field(parsdata0, field_name)
values = [];
for k = 1 : numel(parsdata0)
    if isfield(parsdata0(k), field_name) && ~isempty(parsdata0(k).(field_name))
        values = [values; parsdata0(k).(field_name)];
    end
end
end

function ensdata0 = trim_trial_ensemble(ensdata0, next_row, trial_idx)
if next_row > size(ensdata0.t, 1)
    return
end
ensdata0.t(next_row:end, trial_idx) = 0;
ensdata0.n_agg(next_row:end, trial_idx) = 0;
col = 2 * (trial_idx - 1) + (1:2);
ensdata0.tau(next_row:end, col) = 0;
ensdata0.kn_kin(next_row:end, col) = 0;
ensdata0.kn_diff(next_row:end, col) = 0;
ensdata0.dpp(next_row:end, col) = 0;
ensdata0.sigmapp(next_row:end, col) = 0;
ensdata0.da(next_row:end, col) = 0;
ensdata0.dm(next_row:end, col) = 0;
end

function ensdata0 = trim_empty_ensemble_rows(ensdata0)
keep_rows = any(ensdata0.n_agg > 0, 2);
if ~any(keep_rows)
    return
end
ensdata0.t = ensdata0.t(keep_rows, :);
ensdata0.n_agg = ensdata0.n_agg(keep_rows, :);
ensdata0.tau = ensdata0.tau(keep_rows, :);
ensdata0.kn_kin = ensdata0.kn_kin(keep_rows, :);
ensdata0.kn_diff = ensdata0.kn_diff(keep_rows, :);
ensdata0.dpp = ensdata0.dpp(keep_rows, :);
ensdata0.sigmapp = ensdata0.sigmapp(keep_rows, :);
ensdata0.da = ensdata0.da(keep_rows, :);
ensdata0.dm = ensdata0.dm(keep_rows, :);
end
