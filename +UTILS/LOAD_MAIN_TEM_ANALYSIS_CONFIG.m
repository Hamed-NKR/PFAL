function cfg = LOAD_MAIN_TEM_ANALYSIS_CONFIG(config_path)
%LOAD_MAIN_TEM_ANALYSIS_CONFIG Load and validate the TEM analysis config.
%   CFG = UTILS.LOAD_MAIN_TEM_ANALYSIS_CONFIG reads
%   PFAL_MAIN_TEM_ANALYSIS_CONFIG when set; otherwise it reads
%   config/main_tem_analysis/main_tem_analysis_config.local.json.
%
%   CFG = UTILS.LOAD_MAIN_TEM_ANALYSIS_CONFIG(CONFIG_PATH) reads a
%   specific config file instead of the environment/default-selected file.

if nargin < 1 || isempty(config_path)
    config_path = getenv('PFAL_MAIN_TEM_ANALYSIS_CONFIG');
    if isempty(config_path)
        config_path = fullfile(repo_root(), 'config', 'main_tem_analysis', ...
            'main_tem_analysis_config.local.json');
    end
end

config_path = char(config_path);

if ~isfile(config_path)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:MissingConfig', ...
        ['TEM analysis config not found: %s\n' ...
        'Create it from config/main_tem_analysis/main_tem_analysis_config.example.json, ' ...
        'or set PFAL_MAIN_TEM_ANALYSIS_CONFIG to an explicit config path.'], ...
        config_path);
end

try
    cfg = jsondecode(fileread(config_path));
catch err
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidConfig', ...
        'Could not parse TEM analysis config "%s": %s', ...
        config_path, err.message);
end

config_dir = fileparts(config_path);
cfg.config_file = config_path;

cfg.analysis = require_struct_field(cfg, 'analysis', 'analysis section');
cfg.analysis.coverage_threshold = require_positive_scalar_field( ...
    cfg.analysis, 'coverage_threshold', 'analysis.coverage_threshold');
cfg.analysis.logistic_bandwidth = require_positive_scalar_field( ...
    cfg.analysis, 'logistic_bandwidth', 'analysis.logistic_bandwidth');
cfg.analysis.primary_area_start_row = require_positive_integer_field( ...
    cfg.analysis, 'primary_area_start_row', 'analysis.primary_area_start_row');
cfg.analysis.random_resample_count = require_positive_integer_field( ...
    cfg.analysis, 'random_resample_count', 'analysis.random_resample_count');
cfg.analysis.random_seed = optional_nonnegative_integer_field( ...
    cfg.analysis, 'random_seed');

cfg.outputs = require_struct_field(cfg, 'outputs', 'outputs section');
cfg.outputs.data_root = resolve_config_path(require_text_field( ...
    cfg.outputs, 'data_root', 'outputs.data_root'), config_dir);
cfg.outputs.results_root = resolve_config_path(require_text_field( ...
    cfg.outputs, 'results_root', 'outputs.results_root'), config_dir);
cfg.outputs.mat_file = require_text_field(cfg.outputs, ...
    'mat_file', 'outputs.mat_file');
cfg.outputs.model_inputs_csv = require_text_field(cfg.outputs, ...
    'model_inputs_csv', 'outputs.model_inputs_csv');
cfg.outputs.summary_csv = require_text_field(cfg.outputs, ...
    'summary_csv', 'outputs.summary_csv');

cfg.plots = require_struct_field(cfg, 'plots', 'plots section');
cfg.plots.enabled = require_logical_field(cfg.plots, ...
    'enabled', 'plots.enabled');
cfg.plots.visible = optional_text_field(cfg.plots, 'visible', 'on');
cfg.plots.export = require_logical_field(cfg.plots, ...
    'export', 'plots.export');
cfg.plots.format = optional_text_field(cfg.plots, 'format', 'png');
cfg.plots.resolution = require_positive_integer_field(cfg.plots, ...
    'resolution', 'plots.resolution');
cfg.plots.save_figures = optional_logical_field(cfg.plots, ...
    'save_figures', false);

if ~isfield(cfg, 'entries') || isempty(cfg.entries) || ~isstruct(cfg.entries)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:MissingEntries', ...
        'The TEM analysis config must contain a non-empty entries array.');
end

entries_in = cfg.entries;
entries_out = cell(numel(entries_in), 1);
for i = 1:numel(entries_in)
    entries_out{i} = validate_entry(entries_in(i), i, config_dir);
end
cfg.entries = vertcat(entries_out{:});

end

function entry = validate_entry(entry, index, config_dir)
label = sprintf('entries(%d)', index);

entry.id = require_text_field(entry, 'id', [label '.id']);
entry.label = require_text_field(entry, 'label', [label '.label']);
entry.entry_type = optional_text_field(entry, 'entry_type', '');
entry.entry_dir = resolve_config_path(require_text_field(entry, ...
    'entry_dir', [label '.entry_dir']), config_dir);

entry.aggregate = require_struct_field(entry, 'aggregate', ...
    [label '.aggregate section']);
entry.aggregate.file = resolve_entry_path(require_text_field( ...
    entry.aggregate, 'file', [label '.aggregate.file']), entry.entry_dir);
entry.aggregate.variable = require_text_field(entry.aggregate, ...
    'variable', [label '.aggregate.variable']);

entry.primary_particles = require_struct_field(entry, 'primary_particles', ...
    [label '.primary_particles section']);
entry.primary_particles.folder = resolve_entry_path(require_text_field( ...
    entry.primary_particles, 'folder', [label '.primary_particles.folder']), ...
    entry.entry_dir);
entry.primary_particles.file_pattern = require_text_field( ...
    entry.primary_particles, 'file_pattern', ...
    [label '.primary_particles.file_pattern']);
entry.primary_particles.area_column = require_text_field( ...
    entry.primary_particles, 'area_column', ...
    [label '.primary_particles.area_column']);

if ~isfield(entry, 'aggregate_ids') || isempty(entry.aggregate_ids)
    entry.aggregate_ids = [];
else
    entry.aggregate_ids = parse_aggregate_ids(entry.aggregate_ids, label);
end

entry.plot = optional_struct_field(entry, 'plot');
entry.plot.color = optional_text_field(entry.plot, 'color', '');
entry.plot.marker = optional_text_field(entry.plot, 'marker', '');

entry.include_in_model_inputs = optional_logical_field(entry, ...
    'include_in_model_inputs', true);
entry.include_in_publication_plots = optional_logical_field(entry, ...
    'include_in_publication_plots', true);
entry.include_in_hybridity_scatter = optional_logical_field(entry, ...
    'include_in_hybridity_scatter', true);

if isfield(entry, 'model_filter') && ~isempty(entry.model_filter)
    entry.model_filter = validate_filter(entry.model_filter, label);
else
    entry.model_filter = struct();
end

end

function filter = validate_filter(filter, label)
if ~isstruct(filter)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidModelFilter', ...
        '%s.model_filter must be a struct when provided.', label);
end
filter.field = require_text_field(filter, 'field', ...
    [label '.model_filter.field']);
filter.operator = require_text_field(filter, 'operator', ...
    [label '.model_filter.operator']);
filter.value = require_scalar_field(filter, 'value', ...
    [label '.model_filter.value']);
end

function ids = parse_aggregate_ids(src, label)
if isnumeric(src)
    ids = double(src(:));
elseif isstruct(src)
    if isfield(src, 'values')
        ids = double(src.values(:));
    else
        start_id = require_positive_integer_field(src, 'start', ...
            [label '.aggregate_ids.start']);
        end_id = require_positive_integer_field(src, 'end', ...
            [label '.aggregate_ids.end']);
        step_id = optional_positive_integer_field(src, 'step', 1);
        if end_id < start_id
            error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidAggregateIdRange', ...
                '%s.aggregate_ids.end must be greater than or equal to start.', ...
                label);
        end
        ids = (start_id:step_id:end_id).';
    end
else
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidAggregateIds', ...
        '%s.aggregate_ids must be a numeric array or range struct.', label);
end

if isempty(ids) || any(~isfinite(ids)) || any(ids <= 0) || any(round(ids) ~= ids)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidAggregateIds', ...
        '%s.aggregate_ids must contain positive integer ids.', label);
end
ids = unique(ids(:), 'stable');
end

function out = require_struct_field(src, field_name, label)
if ~isfield(src, field_name) || ~isstruct(src.(field_name))
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:MissingSection', ...
        'The TEM analysis config is missing the %s.', label);
end
out = src.(field_name);
end

function out = optional_struct_field(src, field_name)
if isfield(src, field_name) && ~isempty(src.(field_name))
    if ~isstruct(src.(field_name))
        error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidStruct', ...
            'The TEM analysis config field %s must be a struct.', field_name);
    end
    out = src.(field_name);
else
    out = struct();
end
end

function out = require_text_field(src, field_name, label)
if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:MissingText', ...
        'The TEM analysis config is missing %s.', label);
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
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:MissingScalar', ...
        'The TEM analysis config is missing %s.', label);
end
out = double(src.(field_name));
if ~isscalar(out) || ~isfinite(out)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidScalar', ...
        'The TEM analysis config field %s must be a finite numeric scalar.', ...
        label);
end
end

function out = require_positive_scalar_field(src, field_name, label)
out = require_scalar_field(src, field_name, label);
if out <= 0
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidPositiveScalar', ...
        'The TEM analysis config field %s must be positive.', label);
end
end

function out = require_positive_integer_field(src, field_name, label)
out = require_positive_scalar_field(src, field_name, label);
if round(out) ~= out
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidInteger', ...
        'The TEM analysis config field %s must be a positive integer.', ...
        label);
end
end

function out = optional_positive_integer_field(src, field_name, default_value)
out = default_value;
if isfield(src, field_name) && ~isempty(src.(field_name))
    out = require_positive_integer_field(src, field_name, field_name);
end
end

function out = optional_nonnegative_integer_field(src, field_name)
out = [];
if ~isfield(src, field_name) || isempty(src.(field_name))
    return
end
out = double(src.(field_name));
if ~isscalar(out) || ~isfinite(out) || out < 0 || round(out) ~= out
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidNonnegativeInteger', ...
        'analysis.%s must be empty or a nonnegative integer scalar.', ...
        field_name);
end
end

function out = require_logical_field(src, field_name, label)
if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:MissingLogical', ...
        'The TEM analysis config is missing %s.', label);
end
out = logical(src.(field_name));
if ~isscalar(out)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidLogical', ...
        'The TEM analysis config field %s must be a scalar logical.', label);
end
end

function out = optional_logical_field(src, field_name, default_value)
out = default_value;
if isfield(src, field_name) && ~isempty(src.(field_name))
    out = logical(src.(field_name));
    if ~isscalar(out)
        error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidLogical', ...
            'The TEM analysis config field %s must be a scalar logical.', ...
            field_name);
    end
end
end

function full_path = resolve_entry_path(path_value, entry_dir)
full_path = char(path_value);
if is_absolute_path(full_path)
    full_path = normalize_path(full_path);
    return
end
full_path = normalize_path(fullfile(entry_dir, full_path));
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
