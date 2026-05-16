function cfg = LOAD_MAIN_LD2_CONFIG(config_path)
%LOAD_MAIN_LD2_CONFIG Load and validate the main_LD2 JSON config.
%   CFG = UTILS.LOAD_MAIN_LD2_CONFIG reads the local JSON config for
%   main_LD2_v3 and returns a validated MATLAB struct.
%
%   CFG = UTILS.LOAD_MAIN_LD2_CONFIG(CONFIG_PATH) reads a specific config
%   file instead of the default local file.

if nargin < 1 || isempty(config_path)
    config_path = fullfile(repo_root(), 'config', 'main_ld2_config.local.json');
end

config_path = char(config_path);

if ~isfile(config_path)
    error('PFAL:LOAD_MAIN_LD2_CONFIG:MissingConfig', ...
        ['LD2 config not found: %s\n' ...
        'Create it from config/main_ld2_config.example.json and fill in the local dataset path.'], ...
        config_path);
end

try
    cfg = jsondecode(fileread(config_path));
catch err
    error('PFAL:LOAD_MAIN_LD2_CONFIG:InvalidConfig', ...
        'Could not parse LD2 config "%s": %s', config_path, err.message);
end

config_dir = fileparts(config_path);
cfg.config_file = config_path;

cfg.dataset = require_struct_field(cfg, 'dataset', 'dataset section');
cfg.dataset.id = require_text_field(cfg.dataset, 'id', 'dataset.id');
cfg.dataset.file = require_text_field(cfg.dataset, 'file', 'dataset.file');
cfg.dataset.variable = require_text_field(cfg.dataset, 'variable', 'dataset.variable');
cfg.dataset.url = optional_text_field(cfg.dataset, 'url');
cfg.dataset.doi = optional_text_field(cfg.dataset, 'doi');
cfg.dataset.sha256 = optional_text_field(cfg.dataset, 'sha256');
cfg.dataset.cache_file = optional_text_field(cfg.dataset, 'cache_file');
cfg.dataset.config_dir = config_dir;

cfg.results = require_struct_field(cfg, 'results', 'results section');
cfg.results.root = require_text_field(cfg.results, 'root', 'results.root');
cfg.results.run_label = require_text_field(cfg.results, 'run_label', 'results.run_label');
cfg.results.checkpoint_prefix = require_text_field(cfg.results, ...
    'checkpoint_prefix', 'results.checkpoint_prefix');
cfg.results.final_prefix = require_text_field(cfg.results, ...
    'final_prefix', 'results.final_prefix');
cfg.results.root = resolve_config_path(cfg.results.root, config_dir);

cfg.simulation = require_struct_field(cfg, 'simulation', 'simulation section');
cfg.simulation.k_max = require_positive_integer_field(cfg.simulation, ...
    'k_max', 'simulation.k_max');
cfg.simulation.r_n_agg = require_vector_field(cfg.simulation, ...
    'r_n_agg', 'simulation.r_n_agg');
cfg.simulation.checkpoint_interval = require_positive_integer_field(cfg.simulation, ...
    'checkpoint_interval', 'simulation.checkpoint_interval');

cfg.projection = require_struct_field(cfg, 'projection', 'projection section');
cfg.projection.n_mc = require_positive_integer_field(cfg.projection, ...
    'n_mc', 'projection.n_mc');
cfg.projection.n_ang = require_positive_integer_field(cfg.projection, ...
    'n_ang', 'projection.n_ang');

cfg.transport = require_struct_field(cfg, 'transport', 'transport section');
cfg.transport.user_defined = require_struct_field(cfg.transport, ...
    'user_defined', 'transport.user_defined section');
cfg.transport.fluid_options = require_struct_field(cfg.transport, ...
    'fluid_options', 'transport.fluid_options section');
cfg.transport.mobility_options = require_struct_field(cfg.transport, ...
    'mobility_options', 'transport.mobility_options section');
cfg.transport.location_options = require_struct_field(cfg.transport, ...
    'location_options', 'transport.location_options section');
cfg.transport.growth_options = require_struct_field(cfg.transport, ...
    'growth_options', 'transport.growth_options section');

cfg.transport.fluid_options.amb = require_text_field(cfg.transport.fluid_options, ...
    'amb', 'transport.fluid_options.amb');
cfg.transport.mobility_options.c_dt = require_scalar_field( ...
    cfg.transport.mobility_options, 'c_dt', 'transport.mobility_options.c_dt');
cfg.transport.mobility_options.mtd = require_text_field( ...
    cfg.transport.mobility_options, 'mtd', 'transport.mobility_options.mtd');
cfg.transport.location_options.vf = require_text_field( ...
    cfg.transport.location_options, 'vf', 'transport.location_options.vf');
cfg.transport.growth_options.indupdate = require_text_field( ...
    cfg.transport.growth_options, 'indupdate', 'transport.growth_options.indupdate');

end

function out = require_struct_field(src, field_name, label)
%REQUIRE_STRUCT_FIELD Read and validate a nested struct field.

if ~isfield(src, field_name) || ~isstruct(src.(field_name))
    error('PFAL:LOAD_MAIN_LD2_CONFIG:MissingSection', ...
        'The LD2 config is missing the %s.', label);
end

out = src.(field_name);

end

function out = require_text_field(src, field_name, label)
%REQUIRE_TEXT_FIELD Read and validate a non-empty text field.

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_LD2_CONFIG:MissingText', ...
        'The LD2 config is missing %s.', label);
end

out = char(src.(field_name));

end

function out = optional_text_field(src, field_name)
%OPTIONAL_TEXT_FIELD Read optional text config fields.

out = '';
if isfield(src, field_name) && ~isempty(src.(field_name))
    out = char(src.(field_name));
end

end

function out = require_scalar_field(src, field_name, label)
%REQUIRE_SCALAR_FIELD Read and validate a numeric scalar field.

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_LD2_CONFIG:MissingScalar', ...
        'The LD2 config is missing %s.', label);
end

out = double(src.(field_name));
if ~isscalar(out) || ~isfinite(out)
    error('PFAL:LOAD_MAIN_LD2_CONFIG:InvalidScalar', ...
        'The LD2 config field %s must be a finite numeric scalar.', label);
end

end

function out = require_positive_integer_field(src, field_name, label)
%REQUIRE_POSITIVE_INTEGER_FIELD Read and validate a positive integer field.

out = require_scalar_field(src, field_name, label);
if out <= 0 || round(out) ~= out
    error('PFAL:LOAD_MAIN_LD2_CONFIG:InvalidInteger', ...
        'The LD2 config field %s must be a positive integer.', label);
end

end

function out = require_vector_field(src, field_name, label)
%REQUIRE_VECTOR_FIELD Read and validate a finite numeric vector field.

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_LD2_CONFIG:MissingVector', ...
        'The LD2 config is missing %s.', label);
end

out = double(src.(field_name));
out = reshape(out, 1, []);
if isempty(out) || any(~isfinite(out))
    error('PFAL:LOAD_MAIN_LD2_CONFIG:InvalidVector', ...
        'The LD2 config field %s must be a finite numeric vector.', label);
end

end

function full_path = resolve_config_path(path_value, config_dir)
%RESOLVE_CONFIG_PATH Expand relative paths from the config directory.

full_path = char(path_value);
if is_absolute_path(full_path)
    full_path = normalize_path(full_path);
    return
end

full_path = normalize_path(fullfile(config_dir, full_path));

end

function tf = is_absolute_path(path_value)
%IS_ABSOLUTE_PATH Detect Windows drive, UNC, and POSIX absolute paths.

tf = ~isempty(regexp(path_value, '^[A-Za-z]:[\\/]', 'once')) || ...
    strncmp(path_value, '\\', 2) || ...
    strncmp(path_value, '//', 2) || ...
    strncmp(path_value, '/', 1);

end

function path_value = normalize_path(path_value)
%NORMALIZE_PATH Collapse relative path segments when MATLAB has Java.

try
    path_value = char(java.io.File(path_value).getCanonicalPath());
catch
    path_value = char(path_value);
end

end

function root_dir = repo_root()
%REPO_ROOT Resolve the repository root from the current package location.

utils_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(utils_dir);

end
