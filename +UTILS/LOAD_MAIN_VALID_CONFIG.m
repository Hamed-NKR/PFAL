function cfg = LOAD_MAIN_VALID_CONFIG(config_path)
%LOAD_MAIN_VALID_CONFIG Load and validate the publication validation config.
%   The PFAL_MAIN_VALID_CONFIG environment variable selects a profile when
%   CONFIG_PATH is omitted. Relative paths are resolved from the config file,
%   which keeps machine-specific addresses out of tracked MATLAB source.

if nargin < 1 || isempty(config_path)
    config_path = getenv('PFAL_MAIN_VALID_CONFIG');
    if isempty(config_path)
        config_path = fullfile(repo_root(), 'config', 'main_valid', ...
            'main_valid_config.local.json');
    end
end
config_path = char(config_path);

if ~isfile(config_path)
    error('PFAL:LOAD_MAIN_VALID_CONFIG:MissingConfig', ...
        ['Validation config not found: %s\nCreate it from ' ...
         'config/main_valid/main_valid_config.example.json or set ' ...
         'PFAL_MAIN_VALID_CONFIG.'], config_path);
end

try
    cfg = jsondecode(fileread(config_path));
catch err
    error('PFAL:LOAD_MAIN_VALID_CONFIG:InvalidJSON', ...
        'Could not parse validation config "%s": %s', config_path, err.message);
end

config_dir = fileparts(config_path);
cfg.config_file = normalize_path(config_path);

cfg.sources = require_struct(cfg, 'sources', 'sources');
source_names = {'ld2', 'tem', 'effective_density'};
for i = 1:numel(source_names)
    name = source_names{i};
    cfg.sources.(name) = validate_source( ...
        require_struct(cfg.sources, name, ['sources.', name]), ...
        ['sources.', name], config_dir);
end

if ~isfield(cfg, 'conditions') || isempty(cfg.conditions) || ...
        ~isstruct(cfg.conditions)
    error('PFAL:LOAD_MAIN_VALID_CONFIG:MissingConditions', ...
        'The validation config must define a non-empty conditions array.');
end
condition_cells = cell(numel(cfg.conditions), 1);
for i = 1:numel(cfg.conditions)
    condition_cells{i} = validate_condition(cfg.conditions(i), i);
end
cfg.conditions = vertcat(condition_cells{:});
condition_ids = string({cfg.conditions.id});
if numel(unique(condition_ids)) ~= numel(condition_ids)
    error('PFAL:LOAD_MAIN_VALID_CONFIG:DuplicateConditionID', ...
        'Each configured condition must have a unique id.');
end

cfg.physics = require_struct(cfg, 'physics', 'physics');
cfg.physics.material_density_kg_m3 = positive_scalar(cfg.physics, ...
    'material_density_kg_m3', 'physics.material_density_kg_m3');

cfg.references = require_struct(cfg, 'references', 'references');
cfg.references.dpp_vs_da = validate_reference( ...
    require_struct(cfg.references, 'dpp_vs_da', 'references.dpp_vs_da'), ...
    'references.dpp_vs_da', 'exponent', 'value_at_100_nm');
cfg.references.rho_eff_vs_dm = validate_reference( ...
    require_struct(cfg.references, 'rho_eff_vs_dm', ...
    'references.rho_eff_vs_dm'), 'references.rho_eff_vs_dm', ...
    'mass_mobility_exponent', 'value_at_100_nm');

cfg.fits = require_struct(cfg, 'fits', 'fits');
cfg.fits.dpp_vs_da = validate_fit(require_struct(cfg.fits, ...
    'dpp_vs_da', 'fits.dpp_vs_da'), 'fits.dpp_vs_da');
cfg.fits.rho_eff_vs_dm = validate_fit(require_struct(cfg.fits, ...
    'rho_eff_vs_dm', 'fits.rho_eff_vs_dm'), 'fits.rho_eff_vs_dm');

cfg.figures = validate_figures(require_struct(cfg, 'figures', 'figures'));
cfg.outputs = require_struct(cfg, 'outputs', 'outputs');
cfg.outputs.root = resolve_config_path(require_text(cfg.outputs, ...
    'root', 'outputs.root'), config_dir);
cfg.outputs.export = logical_scalar(cfg.outputs, 'export', 'outputs.export');
cfg.outputs.pdf = logical_scalar(cfg.outputs, 'pdf', 'outputs.pdf');
cfg.outputs.png = logical_scalar(cfg.outputs, 'png', 'outputs.png');
cfg.outputs.png_resolution = positive_integer(cfg.outputs, ...
    'png_resolution', 'outputs.png_resolution');
cfg.outputs.save_figures = logical_scalar(cfg.outputs, ...
    'save_figures', 'outputs.save_figures');

end

function source = validate_source(source, label, config_dir)
%VALIDATE_SOURCE Normalize source metadata and config-relative paths.

source.id = require_text(source, 'id', [label, '.id']);
source.file = resolve_config_path(require_text(source, 'file', ...
    [label, '.file']), config_dir);
source.url = optional_text(source, 'url', '');
source.doi = optional_text(source, 'doi', '');
source.sha256 = optional_text(source, 'sha256', '');
source.cache_file = optional_text(source, 'cache_file', '');
if ~isempty(source.cache_file)
    source.cache_file = resolve_config_path(source.cache_file, config_dir);
end
if ~isfield(source, 'variables') || isempty(source.variables)
    error('PFAL:LOAD_MAIN_VALID_CONFIG:MissingSourceVariables', ...
        '%s.variables must list the required MAT variables.', label);
end

end

function condition = validate_condition(condition, index)
%VALIDATE_CONDITION Validate named experimental/simulation mappings.

label = sprintf('conditions(%d)', index);
condition.id = require_text(condition, 'id', [label, '.id']);
condition.label = require_text(condition, 'label', [label, '.label']);
condition.enabled = logical_scalar(condition, 'enabled', [label, '.enabled']);
condition.include = require_struct(condition, 'include', [label, '.include']);
condition.include.dpp_vs_da = logical_scalar(condition.include, ...
    'dpp_vs_da', [label, '.include.dpp_vs_da']);
condition.include.rho_eff_vs_dm = logical_scalar(condition.include, ...
    'rho_eff_vs_dm', [label, '.include.rho_eff_vs_dm']);
condition.tem_entry_id = optional_text(condition, 'tem_entry_id', '');
condition.effective_density_condition_id = optional_text(condition, ...
    'effective_density_condition_id', '');
condition.ld2_fractions = optional_numeric_vector(condition, ...
    'ld2_fractions', []);
condition.style = require_struct(condition, 'style', [label, '.style']);
condition.style.color = require_text(condition.style, 'color', ...
    [label, '.style.color']);
condition.style.simulation_marker = require_text(condition.style, ...
    'simulation_marker', [label, '.style.simulation_marker']);
condition.style.experimental_marker = require_text(condition.style, ...
    'experimental_marker', [label, '.style.experimental_marker']);
condition.style.simulation_marker_size = positive_scalar(condition.style, ...
    'simulation_marker_size', [label, '.style.simulation_marker_size']);
condition.style.experimental_marker_size = positive_scalar(condition.style, ...
    'experimental_marker_size', [label, '.style.experimental_marker_size']);
condition.style.fit_line_width = positive_scalar(condition.style, ...
    'fit_line_width', [label, '.style.fit_line_width']);

if condition.enabled && ...
        (condition.include.dpp_vs_da || condition.include.rho_eff_vs_dm) && ...
        isempty(condition.ld2_fractions)
    error('PFAL:LOAD_MAIN_VALID_CONFIG:MissingLD2Mapping', ...
        '%s must provide ld2_fractions when an enabled panel uses simulation data.', ...
        label);
end
if condition.enabled && condition.include.dpp_vs_da && ...
        isempty(condition.tem_entry_id)
    error('PFAL:LOAD_MAIN_VALID_CONFIG:MissingTEMMapping', ...
        '%s must provide tem_entry_id for the dpp_vs_da panel.', label);
end
if condition.enabled && condition.include.rho_eff_vs_dm && ...
        isempty(condition.effective_density_condition_id)
    error('PFAL:LOAD_MAIN_VALID_CONFIG:MissingDensityMapping', ...
        '%s must provide effective_density_condition_id for rho_eff_vs_dm.', label);
end

end

function reference = validate_reference(reference, label, exponent_field, value_field)
%VALIDATE_REFERENCE Validate a configurable literature correlation.

reference.enabled = logical_scalar(reference, 'enabled', [label, '.enabled']);
reference.label = require_text(reference, 'label', [label, '.label']);
reference.(exponent_field) = finite_scalar(reference, exponent_field, ...
    [label, '.', exponent_field]);
reference.(value_field) = positive_scalar(reference, value_field, ...
    [label, '.', value_field]);
reference.color = require_text(reference, 'color', [label, '.color']);
reference.line_style = require_text(reference, 'line_style', ...
    [label, '.line_style']);
reference.line_width = positive_scalar(reference, 'line_width', ...
    [label, '.line_width']);

end

function fit_cfg = validate_fit(fit_cfg, label)
%VALIDATE_FIT Validate independent Bayesian controls for one panel.

fit_cfg.enabled = logical_scalar(fit_cfg, 'enabled', [label, '.enabled']);
fit_cfg.degree = positive_integer(fit_cfg, 'degree', [label, '.degree']);
fit_cfg.resolution = positive_integer(fit_cfg, 'resolution', ...
    [label, '.resolution']);
if fit_cfg.resolution < 2
    error('PFAL:LOAD_MAIN_VALID_CONFIG:InvalidResolution', ...
        '%s.resolution must be at least 2.', label);
end
fit_cfg.posterior_samples = positive_integer(fit_cfg, ...
    'posterior_samples', [label, '.posterior_samples']);
if fit_cfg.posterior_samples < 100
    error('PFAL:LOAD_MAIN_VALID_CONFIG:InvalidPosteriorSamples', ...
        '%s.posterior_samples must be at least 100.', label);
end
fit_cfg.rng_seed = nonnegative_integer(fit_cfg, 'rng_seed', ...
    [label, '.rng_seed']);
fit_cfg.credible_level = finite_scalar(fit_cfg, 'credible_level', ...
    [label, '.credible_level']);
if fit_cfg.credible_level <= 0 || fit_cfg.credible_level >= 1
    error('PFAL:LOAD_MAIN_VALID_CONFIG:InvalidCredibleLevel', ...
        '%s.credible_level must be strictly between 0 and 1.', label);
end
fit_cfg.show_band = logical_scalar(fit_cfg, 'show_band', ...
    [label, '.show_band']);
fit_cfg.x_transform = validate_transform(require_text(fit_cfg, ...
    'x_transform', [label, '.x_transform']), [label, '.x_transform']);
fit_cfg.y_transform = validate_transform(require_text(fit_cfg, ...
    'y_transform', [label, '.y_transform']), [label, '.y_transform']);
fit_cfg.slope_offset = finite_scalar(fit_cfg, 'slope_offset', ...
    [label, '.slope_offset']);
fit_cfg.prior_mu = optional_numeric_vector(fit_cfg, 'prior_mu', []);
fit_cfg.prior_v_scale = positive_scalar(fit_cfg, 'prior_v_scale', ...
    [label, '.prior_v_scale']);
fit_cfg.prior_a = positive_scalar(fit_cfg, 'prior_a', [label, '.prior_a']);
fit_cfg.prior_b = positive_scalar(fit_cfg, 'prior_b', [label, '.prior_b']);

end

function figures = validate_figures(figures)
%VALIDATE_FIGURES Validate shared typography and panel-specific appearance.

figures.visible = lower(require_text(figures, 'visible', 'figures.visible'));
if ~ismember(figures.visible, {'on', 'off'})
    error('PFAL:LOAD_MAIN_VALID_CONFIG:InvalidVisibility', ...
        'figures.visible must be "on" or "off".');
end
figures.show_raw_simulation = logical_scalar(figures, ...
    'show_raw_simulation', 'figures.show_raw_simulation');
figures.show_tem_error_bars = logical_scalar(figures, ...
    'show_tem_error_bars', 'figures.show_tem_error_bars');
figures.simulation_marker_alpha = unit_interval(figures, ...
    'simulation_marker_alpha', 'figures.simulation_marker_alpha');
figures.band_alpha = unit_interval(figures, 'band_alpha', ...
    'figures.band_alpha');
figures.font = require_struct(figures, 'font', 'figures.font');
figures.font.family = require_text(figures.font, 'family', ...
    'figures.font.family');
figures.font.fallback = require_text(figures.font, 'fallback', ...
    'figures.font.fallback');
figures.font.axis_size = positive_scalar(figures.font, ...
    'axis_size', 'figures.font.axis_size');
figures.font.label_size = positive_scalar(figures.font, ...
    'label_size', 'figures.font.label_size');
figures.font.legend_size = positive_scalar(figures.font, ...
    'legend_size', 'figures.font.legend_size');
figures.font.weight = require_text(figures.font, 'weight', ...
    'figures.font.weight');
figures.font.interpreter = lower(require_text(figures.font, ...
    'interpreter', 'figures.font.interpreter'));
if ~ismember(figures.font.interpreter, {'tex', 'none', 'latex'})
    error('PFAL:LOAD_MAIN_VALID_CONFIG:InvalidInterpreter', ...
        'figures.font.interpreter must be tex, none, or latex.');
end

panel_names = {'dpp_vs_da', 'rho_eff_vs_dm'};
for i = 1:numel(panel_names)
    name = panel_names{i};
    panel = require_struct(figures, name, ['figures.', name]);
    panel.position = numeric_vector(panel, 'position', 4, ...
        ['figures.', name, '.position']);
    panel.x_limits = optional_limits(panel, 'x_limits', ...
        ['figures.', name, '.x_limits']);
    panel.y_limits = optional_limits(panel, 'y_limits', ...
        ['figures.', name, '.y_limits']);
    figures.(name) = panel;
end

end

function value = validate_transform(value, label)
%VALIDATE_TRANSFORM Restrict fit transforms to FIT_POLY-supported values.

value = lower(value);
if ~ismember(value, {'log10', 'none'})
    error('PFAL:LOAD_MAIN_VALID_CONFIG:InvalidTransform', ...
        '%s must be "log10" or "none".', label);
end

end

function out = require_struct(source, field_name, label)
if ~isfield(source, field_name) || ~isstruct(source.(field_name))
    error('PFAL:LOAD_MAIN_VALID_CONFIG:MissingSection', ...
        'The validation config is missing the %s section.', label);
end
out = source.(field_name);
end

function out = require_text(source, field_name, label)
if ~isfield(source, field_name) || isempty(source.(field_name))
    error('PFAL:LOAD_MAIN_VALID_CONFIG:MissingText', ...
        'The validation config is missing %s.', label);
end
out = char(source.(field_name));
end

function out = optional_text(source, field_name, default_value)
out = default_value;
if isfield(source, field_name) && ~isempty(source.(field_name))
    out = char(source.(field_name));
end
end

function out = finite_scalar(source, field_name, label)
if ~isfield(source, field_name) || isempty(source.(field_name))
    error('PFAL:LOAD_MAIN_VALID_CONFIG:MissingNumber', ...
        'The validation config is missing %s.', label);
end
out = double(source.(field_name));
if ~isscalar(out) || ~isfinite(out)
    error('PFAL:LOAD_MAIN_VALID_CONFIG:InvalidNumber', ...
        '%s must be a finite numeric scalar.', label);
end
end

function out = positive_scalar(source, field_name, label)
out = finite_scalar(source, field_name, label);
if out <= 0
    error('PFAL:LOAD_MAIN_VALID_CONFIG:InvalidPositiveNumber', ...
        '%s must be greater than zero.', label);
end
end

function out = positive_integer(source, field_name, label)
out = positive_scalar(source, field_name, label);
if out ~= round(out)
    error('PFAL:LOAD_MAIN_VALID_CONFIG:InvalidPositiveInteger', ...
        '%s must be a positive integer.', label);
end
end

function out = nonnegative_integer(source, field_name, label)
out = finite_scalar(source, field_name, label);
if out < 0 || out ~= round(out)
    error('PFAL:LOAD_MAIN_VALID_CONFIG:InvalidNonnegativeInteger', ...
        '%s must be a nonnegative integer.', label);
end
end

function out = logical_scalar(source, field_name, label)
if ~isfield(source, field_name) || ~isscalar(source.(field_name)) || ...
        ~islogical(source.(field_name))
    error('PFAL:LOAD_MAIN_VALID_CONFIG:InvalidLogical', ...
        '%s must be true or false.', label);
end
out = source.(field_name);
end

function out = unit_interval(source, field_name, label)
out = finite_scalar(source, field_name, label);
if out < 0 || out > 1
    error('PFAL:LOAD_MAIN_VALID_CONFIG:InvalidUnitInterval', ...
        '%s must be between 0 and 1.', label);
end
end

function out = optional_numeric_vector(source, field_name, default_value)
out = default_value;
if isfield(source, field_name) && ~isempty(source.(field_name))
    out = double(source.(field_name));
    out = reshape(out, 1, []);
    if any(~isfinite(out))
        error('PFAL:LOAD_MAIN_VALID_CONFIG:InvalidVector', ...
            '%s must contain only finite numeric values.', field_name);
    end
end
end

function out = numeric_vector(source, field_name, expected_count, label)
out = optional_numeric_vector(source, field_name, []);
if numel(out) ~= expected_count
    error('PFAL:LOAD_MAIN_VALID_CONFIG:InvalidVectorLength', ...
        '%s must contain %d numeric values.', label, expected_count);
end
end

function out = optional_limits(source, field_name, label)
out = optional_numeric_vector(source, field_name, []);
if ~isempty(out) && (numel(out) ~= 2 || any(out <= 0) || out(2) <= out(1))
    error('PFAL:LOAD_MAIN_VALID_CONFIG:InvalidLimits', ...
        '%s must be empty or contain two increasing positive values.', label);
end
end

function full_path = resolve_config_path(path_value, config_dir)
if is_absolute_path(path_value)
    full_path = normalize_path(path_value);
else
    full_path = normalize_path(fullfile(config_dir, path_value));
end
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
