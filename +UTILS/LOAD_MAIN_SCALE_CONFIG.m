function cfg = LOAD_MAIN_SCALE_CONFIG(config_path)
%LOAD_MAIN_SCALE_CONFIG Load and validate the main_scale JSON config.
%   CFG = UTILS.LOAD_MAIN_SCALE_CONFIG reads the local JSON config for
%   main_scale_v2 and returns a validated MATLAB struct.
%
%   CFG = UTILS.LOAD_MAIN_SCALE_CONFIG(CONFIG_PATH) reads a specific
%   config file instead of the default local file.

if nargin < 1 || isempty(config_path)
    config_path = fullfile(repo_root(), 'config', 'main_scale', ...
        'main_scale_config.local.json');
end

config_path = char(config_path);

if ~isfile(config_path)
    error('PFAL:LOAD_MAIN_SCALE_CONFIG:MissingConfig', ...
        ['Scale config not found: %s\n' ...
        'Create it from config/main_scale/main_scale_config.example.json and fill in the local dataset path.'], ...
        config_path);
end

try
    cfg = jsondecode(fileread(config_path));
catch err
    error('PFAL:LOAD_MAIN_SCALE_CONFIG:InvalidConfig', ...
        'Could not parse scale config "%s": %s', config_path, err.message);
end

cfg.config_file = config_path;

cfg.dataset = require_struct_field(cfg, 'dataset', 'dataset section');
cfg.dataset.id = require_text_field(cfg.dataset, 'id', 'dataset.id');
cfg.dataset.file = require_text_field(cfg.dataset, 'file', 'dataset.file');
cfg.dataset.variable = require_text_field(cfg.dataset, 'variable', 'dataset.variable');
cfg.dataset.config_dir = fileparts(config_path);

cfg.correlation = require_struct_field(cfg, 'correlation', 'correlation section');

cfg.correlation.universal = require_struct_field(cfg.correlation, ...
    'universal', 'correlation.universal section');
cfg.correlation.universal.D_TEM = require_scalar_field(cfg.correlation.universal, ...
    'D_TEM', 'correlation.universal.D_TEM');
cfg.correlation.universal.dpp100 = require_scalar_field(cfg.correlation.universal, ...
    'dpp100', 'correlation.universal.dpp100');
cfg.correlation.universal.da_lim_uc = require_vector_field(cfg.correlation.universal, ...
    'da_lim_uc', 2, 'correlation.universal.da_lim_uc');
cfg.correlation.universal.n_da_uc = require_positive_integer_field( ...
    cfg.correlation.universal, 'n_da_uc', 'correlation.universal.n_da_uc');

cfg.correlation.brasil = require_struct_field(cfg.correlation, ...
    'brasil', 'correlation.brasil section');
cfg.correlation.brasil.alpha_a = require_scalar_field(cfg.correlation.brasil, ...
    'alpha_a', 'correlation.brasil.alpha_a');
cfg.correlation.brasil.k_a = require_scalar_field(cfg.correlation.brasil, ...
    'k_a', 'correlation.brasil.k_a');
cfg.correlation.brasil.npp_lim_bc = require_vector_field(cfg.correlation.brasil, ...
    'npp_lim_bc', 2, 'correlation.brasil.npp_lim_bc');
cfg.correlation.brasil.n_npp_uc = require_positive_integer_field( ...
    cfg.correlation.brasil, 'n_npp_uc', 'correlation.brasil.n_npp_uc');

cfg.distribution = require_struct_field(cfg, 'distribution', 'distribution section');
cfg.distribution.gm_da = require_scalar_field(cfg.distribution, ...
    'gm_da', 'distribution.gm_da');
cfg.distribution.gsd_da = require_scalar_field(cfg.distribution, ...
    'gsd_da', 'distribution.gsd_da');
cfg.distribution.gm_dpp = require_scalar_field(cfg.distribution, ...
    'gm_dpp', 'distribution.gm_dpp');
cfg.distribution.gsd_dpp = require_scalar_field(cfg.distribution, ...
    'gsd_dpp', 'distribution.gsd_dpp');

cfg.sampling = require_struct_field(cfg, 'sampling', 'sampling section');
cfg.sampling.cn_scale = require_scalar_field(cfg.sampling, ...
    'cn_scale', 'sampling.cn_scale');
cfg.sampling.n_bin_filter = require_positive_integer_field(cfg.sampling, ...
    'n_bin_filter', 'sampling.n_bin_filter');

cfg.projection = require_struct_field(cfg, 'projection', 'projection section');
cfg.projection.n_mc = require_positive_integer_field(cfg.projection, ...
    'n_mc', 'projection.n_mc');
cfg.projection.n_ang = require_positive_integer_field(cfg.projection, ...
    'n_ang', 'projection.n_ang');

cfg.options = require_struct_field(cfg, 'options', 'options section');

end

function out = require_struct_field(src, field_name, label)
%REQUIRE_STRUCT_FIELD Read and validate a nested struct field.

if ~isfield(src, field_name) || ~isstruct(src.(field_name))
    error('PFAL:LOAD_MAIN_SCALE_CONFIG:MissingSection', ...
        'The scale config is missing the %s.', label);
end

out = src.(field_name);

end

function out = require_text_field(src, field_name, label)
%REQUIRE_TEXT_FIELD Read and validate a non-empty text field.

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_SCALE_CONFIG:MissingText', ...
        'The scale config is missing %s.', label);
end

out = char(src.(field_name));

end

function out = require_scalar_field(src, field_name, label)
%REQUIRE_SCALAR_FIELD Read and validate a numeric scalar field.

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_SCALE_CONFIG:MissingScalar', ...
        'The scale config is missing %s.', label);
end

out = double(src.(field_name));
if ~isscalar(out) || ~isfinite(out)
    error('PFAL:LOAD_MAIN_SCALE_CONFIG:InvalidScalar', ...
        'The scale config field %s must be a finite numeric scalar.', label);
end

end

function out = require_positive_integer_field(src, field_name, label)
%REQUIRE_POSITIVE_INTEGER_FIELD Read and validate a positive integer field.

out = require_scalar_field(src, field_name, label);
if out <= 0 || round(out) ~= out
    error('PFAL:LOAD_MAIN_SCALE_CONFIG:InvalidInteger', ...
        'The scale config field %s must be a positive integer.', label);
end

end

function out = require_vector_field(src, field_name, expected_len, label)
%REQUIRE_VECTOR_FIELD Read and validate a numeric vector field.

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_SCALE_CONFIG:MissingVector', ...
        'The scale config is missing %s.', label);
end

out = double(src.(field_name));
out = reshape(out, 1, []);
if numel(out) ~= expected_len || any(~isfinite(out))
    error('PFAL:LOAD_MAIN_SCALE_CONFIG:InvalidVector', ...
        'The scale config field %s must be a finite numeric vector of length %d.', ...
        label, expected_len);
end

end

function root_dir = repo_root()
%REPO_ROOT Resolve the repository root from the current package location.

utils_dir = fileparts(mfilename('fullpath')); % .../+UTILS
root_dir = fileparts(utils_dir); % repository root above the package folder

end
