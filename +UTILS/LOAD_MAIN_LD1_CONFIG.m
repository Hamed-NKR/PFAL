function cfg = LOAD_MAIN_LD1_CONFIG(config_path)
%LOAD_MAIN_LD1_CONFIG Load and validate the main_LD1 JSON config.
%   CFG = UTILS.LOAD_MAIN_LD1_CONFIG reads PFAL_MAIN_LD1_CONFIG when set;
%   otherwise it reads config/main_ld1/main_ld1_config.local.json.
%
%   CFG = UTILS.LOAD_MAIN_LD1_CONFIG(CONFIG_PATH) reads a specific config
%   file instead of the environment/default-selected file.

if nargin < 1 || isempty(config_path)
    config_path = getenv('PFAL_MAIN_LD1_CONFIG');
    if isempty(config_path)
        config_path = fullfile(repo_root(), 'config', 'main_ld1', ...
            'main_ld1_config.local.json');
    end
end

config_path = char(config_path);

if ~isfile(config_path)
    error('PFAL:LOAD_MAIN_LD1_CONFIG:MissingConfig', ...
        ['LD1 config not found: %s\n' ...
        'Create it from config/main_ld1/main_ld1_config.example.json, ' ...
        'or set PFAL_MAIN_LD1_CONFIG to an explicit config path.'], ...
        config_path);
end

try
    cfg = jsondecode(fileread(config_path));
catch err
    error('PFAL:LOAD_MAIN_LD1_CONFIG:InvalidConfig', ...
        'Could not parse LD1 config "%s": %s', config_path, err.message);
end

config_dir = fileparts(config_path);
cfg.config_file = config_path;

cfg.results = require_struct_field(cfg, 'results', 'results section');
cfg.results.root = require_text_field(cfg.results, 'root', 'results.root');
cfg.results.library_file = require_text_field(cfg.results, ...
    'library_file', 'results.library_file');
cfg.results.run_label = require_text_field(cfg.results, ...
    'run_label', 'results.run_label');
cfg.results.root = resolve_config_path(cfg.results.root, config_dir);

cfg.simulation = require_struct_field(cfg, 'simulation', 'simulation section');
cfg.simulation.n_temporal = require_positive_integer_field(cfg.simulation, ...
    'n_temporal', 'simulation.n_temporal');
cfg.simulation.n_trial = require_positive_integer_field(cfg.simulation, ...
    'n_trial', 'simulation.n_trial');
cfg.simulation.primary_size_gsd_between_trials = require_positive_scalar_field( ...
    cfg.simulation, 'primary_size_gsd_between_trials', ...
    'simulation.primary_size_gsd_between_trials');
cfg.simulation.npp_min = require_positive_integer_field(cfg.simulation, ...
    'npp_min', 'simulation.npp_min');
cfg.simulation.npp_max = require_positive_integer_field(cfg.simulation, ...
    'npp_max', 'simulation.npp_max');
cfg.simulation.j_max = require_positive_integer_field(cfg.simulation, ...
    'j_max', 'simulation.j_max');
cfg.simulation.rng_seed = optional_rng_seed(cfg.simulation);

if cfg.simulation.n_temporal < 2
    error('PFAL:LOAD_MAIN_LD1_CONFIG:InvalidStorageCount', ...
        'simulation.n_temporal must be at least 2.');
end
if cfg.simulation.npp_max < cfg.simulation.npp_min
    error('PFAL:LOAD_MAIN_LD1_CONFIG:InvalidNppRange', ...
        'simulation.npp_max must be greater than or equal to simulation.npp_min.');
end

cfg.projection = require_struct_field(cfg, 'projection', 'projection section');
cfg.projection.n_mc = require_positive_integer_field(cfg.projection, ...
    'n_mc', 'projection.n_mc');
cfg.projection.n_ang = require_positive_integer_field(cfg.projection, ...
    'n_ang', 'projection.n_ang');
if isfield(cfg.projection, 'options') && isstruct(cfg.projection.options)
    cfg.projection.options = cfg.projection.options;
else
    cfg.projection.options = struct();
end
cfg.projection.options.tbar = optional_text_field(cfg.projection.options, ...
    'tbar', 'off');

cfg.transport = require_struct_field(cfg, 'transport', 'transport section');
cfg.transport.domain = require_struct_field(cfg.transport, ...
    'domain', 'transport.domain section');
cfg.transport.domain.volume_fraction = require_scalar_field( ...
    cfg.transport.domain, 'volume_fraction', 'transport.domain.volume_fraction');
cfg.transport.domain.size = require_vector_field(cfg.transport.domain, ...
    'size', 3, 'transport.domain.size');

cfg.transport.particles = require_struct_field(cfg.transport, ...
    'particles', 'transport.particles section');
cfg.transport.particles.n_par = require_positive_integer_field( ...
    cfg.transport.particles, 'n_par', 'transport.particles.n_par');
cfg.transport.particles.n_pp_mean = require_positive_scalar_field( ...
    cfg.transport.particles, 'n_pp_mean', 'transport.particles.n_pp_mean');
cfg.transport.particles.n_pp_std = require_nonnegative_scalar_field( ...
    cfg.transport.particles, 'n_pp_std', 'transport.particles.n_pp_std');
cfg.transport.particles.d_pp_gm = require_positive_scalar_field( ...
    cfg.transport.particles, 'd_pp_gm', 'transport.particles.d_pp_gm');
cfg.transport.particles.d_pp_gsd_between_aggregates = require_nonnegative_scalar_field( ...
    cfg.transport.particles, 'd_pp_gsd_between_aggregates', ...
    'transport.particles.d_pp_gsd_between_aggregates');
cfg.transport.particles.d_pp_gsd_within_aggregates = require_nonnegative_scalar_field( ...
    cfg.transport.particles, 'd_pp_gsd_within_aggregates', ...
    'transport.particles.d_pp_gsd_within_aggregates');

cfg.transport.fluid = require_struct_field(cfg.transport, ...
    'fluid', 'transport.fluid section');
cfg.transport.fluid.temperature = require_positive_scalar_field( ...
    cfg.transport.fluid, 'temperature', 'transport.fluid.temperature');
cfg.transport.fluid.velocity = require_vector_field(cfg.transport.fluid, ...
    'velocity', 3, 'transport.fluid.velocity');
cfg.transport.fluid.pressure = require_positive_scalar_field( ...
    cfg.transport.fluid, 'pressure', 'transport.fluid.pressure');

cfg.transport.fluid_options = require_struct_field(cfg.transport, ...
    'fluid_options', 'transport.fluid_options section');
cfg.transport.fluid_options.amb = require_text_field( ...
    cfg.transport.fluid_options, 'amb', 'transport.fluid_options.amb');

cfg.transport.mobility_options = require_struct_field(cfg.transport, ...
    'mobility_options', 'transport.mobility_options section');
cfg.transport.mobility_options.mtd = require_text_field( ...
    cfg.transport.mobility_options, 'mtd', 'transport.mobility_options.mtd');
cfg.transport.mobility_options.c_dt = require_positive_scalar_field( ...
    cfg.transport.mobility_options, 'c_dt', 'transport.mobility_options.c_dt');

cfg.transport.location_options = require_struct_field(cfg.transport, ...
    'location_options', 'transport.location_options section');
cfg.transport.location_options.vf = require_text_field( ...
    cfg.transport.location_options, 'vf', 'transport.location_options.vf');
cfg.transport.location_options.vf_type = optional_text_field( ...
    cfg.transport.location_options, 'vf_type', 'pp');

if isfield(cfg.transport, 'growth_options') && isstruct(cfg.transport.growth_options)
    cfg.transport.growth_options = cfg.transport.growth_options;
else
    cfg.transport.growth_options = struct();
end

end

function out = require_struct_field(src, field_name, label)
if ~isfield(src, field_name) || ~isstruct(src.(field_name))
    error('PFAL:LOAD_MAIN_LD1_CONFIG:MissingSection', ...
        'The LD1 config is missing the %s.', label);
end
out = src.(field_name);
end

function out = require_text_field(src, field_name, label)
if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_LD1_CONFIG:MissingText', ...
        'The LD1 config is missing %s.', label);
end
out = char(src.(field_name));
end

function out = optional_text_field(src, field_name, default_value)
out = default_value;
if isfield(src, field_name) && ~isempty(src.(field_name))
    out = char(src.(field_name));
end
end

function out = require_scalar_field(src, field_name, label)
if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_LD1_CONFIG:MissingScalar', ...
        'The LD1 config is missing %s.', label);
end
out = double(src.(field_name));
if ~isscalar(out) || ~isfinite(out)
    error('PFAL:LOAD_MAIN_LD1_CONFIG:InvalidScalar', ...
        'The LD1 config field %s must be a finite numeric scalar.', label);
end
end

function out = require_positive_scalar_field(src, field_name, label)
out = require_scalar_field(src, field_name, label);
if out <= 0
    error('PFAL:LOAD_MAIN_LD1_CONFIG:InvalidPositiveScalar', ...
        'The LD1 config field %s must be positive.', label);
end
end

function out = require_nonnegative_scalar_field(src, field_name, label)
out = require_scalar_field(src, field_name, label);
if out < 0
    error('PFAL:LOAD_MAIN_LD1_CONFIG:InvalidNonnegativeScalar', ...
        'The LD1 config field %s must be nonnegative.', label);
end
end

function out = require_positive_integer_field(src, field_name, label)
out = require_positive_scalar_field(src, field_name, label);
if round(out) ~= out
    error('PFAL:LOAD_MAIN_LD1_CONFIG:InvalidInteger', ...
        'The LD1 config field %s must be a positive integer.', label);
end
end

function out = require_vector_field(src, field_name, expected_len, label)
if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_LD1_CONFIG:MissingVector', ...
        'The LD1 config is missing %s.', label);
end
out = double(src.(field_name));
out = reshape(out, 1, []);
if numel(out) ~= expected_len || any(~isfinite(out))
    error('PFAL:LOAD_MAIN_LD1_CONFIG:InvalidVector', ...
        'The LD1 config field %s must be a finite numeric vector of length %d.', ...
        label, expected_len);
end
end

function seed = optional_rng_seed(src)
seed = [];
if ~isfield(src, 'rng_seed') || isempty(src.rng_seed)
    return
end
seed = double(src.rng_seed);
if ~isscalar(seed) || ~isfinite(seed) || seed < 0 || round(seed) ~= seed
    error('PFAL:LOAD_MAIN_LD1_CONFIG:InvalidRngSeed', ...
        'simulation.rng_seed must be empty or a nonnegative integer scalar.');
end
end

function full_path = resolve_config_path(path_value, config_dir)
full_path = char(path_value);
if is_absolute_path(full_path)
    full_path = normalize_path(full_path);
    return
end
full_path = normalize_path(fullfile(config_dir, full_path));
end

function tf = is_absolute_path(path_value)
tf = ~isempty(regexp(path_value, '^[A-Za-z]:[\\/]', 'once')) || ...
    strncmp(path_value, '\\', 2) || strncmp(path_value, '//', 2) || ...
    strncmp(path_value, '/', 1);
end

function path_value = normalize_path(path_value)
try
    path_value = char(java.io.File(path_value).getCanonicalPath());
catch
    path_value = char(path_value);
end
end

function root_dir = repo_root()
utils_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(utils_dir);
end
