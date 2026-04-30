function cfg = LOAD_MAIN_SCATTER_CONFIG(config_path)
%LOAD_MAIN_SCATTER_CONFIG Load and validate the main_scatter JSON config.
%   CFG = UTILS.LOAD_MAIN_SCATTER_CONFIG reads the local JSON config for
%   main_scatter_v8 and returns a validated MATLAB struct.
%
%   CFG = UTILS.LOAD_MAIN_SCATTER_CONFIG(CONFIG_PATH) reads a specific
%   config file instead of the default local file.

if nargin < 1 || isempty(config_path)
    config_path = fullfile(repo_root(), 'config', 'main_scatter_config.local.json');
end

config_path = char(config_path);

% The config file is the single source of user-tuned defaults for the
% scatter workflow, so execution should stop early if it is missing.
if ~isfile(config_path)
    error('PFAL:LOAD_MAIN_SCATTER_CONFIG:MissingConfig', ...
        ['Scatter config not found: %s\n' ...
        'Create it from config/main_scatter_config.example.json and fill in the local dataset path.'], ...
        config_path);
end

% Parse the JSON file into a MATLAB struct before validating each section.
try
    cfg = jsondecode(fileread(config_path));
catch err
    error('PFAL:LOAD_MAIN_SCATTER_CONFIG:InvalidConfig', ...
        'Could not parse scatter config "%s": %s', config_path, err.message);
end

% Validate the dataset section and attach the config directory so the
% dataset loader can resolve relative paths without re-reading the JSON.
cfg.dataset = require_struct_field(cfg, 'dataset', 'dataset section');
cfg.dataset.id = require_text_field(cfg.dataset, 'id', 'dataset.id');
cfg.dataset.file = require_text_field(cfg.dataset, 'file', 'dataset.file');
cfg.dataset.variable = require_text_field(cfg.dataset, 'variable', 'dataset.variable');
cfg.dataset.config_dir = fileparts(config_path);

% Validate the correlation inputs that define the universal and Brasil
% relationships used throughout the script.
cfg.correlation = require_struct_field(cfg, 'correlation', 'correlation section');

cfg.correlation.universal = require_struct_field(cfg.correlation, ...
    'universal', 'correlation.universal section');
cfg.correlation.universal.D_TEM = require_scalar_field(cfg.correlation.universal, ...
    'D_TEM', 'correlation.universal.D_TEM');
cfg.correlation.universal.dpp100 = require_scalar_field(cfg.correlation.universal, ...
    'dpp100', 'correlation.universal.dpp100');
cfg.correlation.universal.da_lim_uc = require_vector_field(cfg.correlation.universal, ...
    'da_lim_uc', 2, 'correlation.universal.da_lim_uc');
cfg.correlation.universal.n_da_uc = require_scalar_field(cfg.correlation.universal, ...
    'n_da_uc', 'correlation.universal.n_da_uc');

cfg.correlation.brasil = require_struct_field(cfg.correlation, ...
    'brasil', 'correlation.brasil section');
cfg.correlation.brasil.alpha_a = require_scalar_field(cfg.correlation.brasil, ...
    'alpha_a', 'correlation.brasil.alpha_a');
cfg.correlation.brasil.k_a = require_scalar_field(cfg.correlation.brasil, ...
    'k_a', 'correlation.brasil.k_a');
cfg.correlation.brasil.npp_lim_bc = require_vector_field(cfg.correlation.brasil, ...
    'npp_lim_bc', 2, 'correlation.brasil.npp_lim_bc');
cfg.correlation.brasil.n_npp_uc = require_scalar_field(cfg.correlation.brasil, ...
    'n_npp_uc', 'correlation.brasil.n_npp_uc');

% Validate the distribution, sampling, projection, and option groups that
% were previously assigned directly in the script.
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
cfg.sampling.cn_scat = require_scalar_field(cfg.sampling, ...
    'cn_scat', 'sampling.cn_scat');

cfg.projection = require_struct_field(cfg, 'projection', 'projection section');
cfg.projection.n_mc = require_scalar_field(cfg.projection, ...
    'n_mc', 'projection.n_mc');
cfg.projection.n_ang = require_scalar_field(cfg.projection, ...
    'n_ang', 'projection.n_ang');

cfg.options = require_struct_field(cfg, 'options', 'options section');
cfg.options.opt_scale = require_text_field(cfg.options, ...
    'opt_scale', 'options.opt_scale');
cfg.options.opts_nppcor = require_text_field(cfg.options, ...
    'opts_nppcor', 'options.opts_nppcor');

end

function out = require_struct_field(src, field_name, label)
%REQUIRE_STRUCT_FIELD Read and validate a nested struct field.

if ~isfield(src, field_name) || ~isstruct(src.(field_name))
    error('PFAL:LOAD_MAIN_SCATTER_CONFIG:MissingSection', ...
        'The scatter config is missing the %s.', label);
end

out = src.(field_name);

end

function out = require_text_field(src, field_name, label)
%REQUIRE_TEXT_FIELD Read and validate a non-empty text field.

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_SCATTER_CONFIG:MissingText', ...
        'The scatter config is missing %s.', label);
end

out = char(src.(field_name));

end

function out = require_scalar_field(src, field_name, label)
%REQUIRE_SCALAR_FIELD Read and validate a numeric scalar field.

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_SCATTER_CONFIG:MissingScalar', ...
        'The scatter config is missing %s.', label);
end

out = double(src.(field_name));
if ~isscalar(out) || ~isfinite(out)
    error('PFAL:LOAD_MAIN_SCATTER_CONFIG:InvalidScalar', ...
        'The scatter config field %s must be a finite numeric scalar.', label);
end

end

function out = require_vector_field(src, field_name, expected_len, label)
%REQUIRE_VECTOR_FIELD Read and validate a numeric vector field.

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_SCATTER_CONFIG:MissingVector', ...
        'The scatter config is missing %s.', label);
end

out = double(src.(field_name));
out = reshape(out, 1, []);
if numel(out) ~= expected_len || any(~isfinite(out))
    error('PFAL:LOAD_MAIN_SCATTER_CONFIG:InvalidVector', ...
        'The scatter config field %s must be a finite numeric vector of length %d.', ...
        label, expected_len);
end

end

function root_dir = repo_root()
%REPO_ROOT Resolve the repository root from the current package location.

utils_dir = fileparts(mfilename('fullpath')); % .../+UTILS
root_dir = fileparts(utils_dir); % repository root above the package folder

end
