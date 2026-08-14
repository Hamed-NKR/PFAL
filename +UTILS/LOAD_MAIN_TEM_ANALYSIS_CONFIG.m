function cfg = LOAD_MAIN_TEM_ANALYSIS_CONFIG(config_path)
%LOAD_MAIN_TEM_ANALYSIS_CONFIG Load and validate the TEM analysis config.
%   CFG = UTILS.LOAD_MAIN_TEM_ANALYSIS_CONFIG reads
%   PFAL_MAIN_TEM_ANALYSIS_CONFIG when set; otherwise it reads
%   config/main_tem_analysis/main_tem_analysis_config.local.json.
%
%   Canonical configs define every publication choice explicitly. Legacy
%   configs are upgraded in memory so older custom profiles remain usable.

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

if ~isfield(cfg, 'entries') || isempty(cfg.entries) || ...
        ~(isstruct(cfg.entries) || iscell(cfg.entries))
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:MissingEntries', ...
        'The TEM analysis config must contain a non-empty entries array.');
end
entries_in = cfg.entries;
entries_out = cell(numel(entries_in), 1);
for i = 1:numel(entries_in)
    if iscell(entries_in)
        entry_in = entries_in{i};
    else
        entry_in = entries_in(i);
    end
    if ~isstruct(entry_in) || ~isscalar(entry_in)
        error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidEntry', ...
            'entries(%d) must be a JSON object.', i);
    end
    entries_out{i} = validate_entry(entry_in, i, config_dir);
end
cfg.entries = vertcat(entries_out{:});
entry_ids = {cfg.entries.id};
if numel(unique(entry_ids)) ~= numel(entry_ids)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:DuplicateEntryId', ...
        'Every entries.id value must be unique.');
end
if ~any([cfg.entries.enabled])
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:NoEnabledEntries', ...
        'At least one TEM entry must be enabled.');
end

cfg.plots = validate_plots(cfg, entry_ids);
cfg.entries = validate_stack_palette_lengths(cfg.entries, ...
    cfg.plots.frequency_bins);

end

function entry = validate_entry(entry, index, config_dir)

% Resolve dataset paths relative to the selected profile and validate the
% complete condition-level scientific and presentation definition.

label = sprintf('entries(%d)', index);
entry.id = require_text_field(entry, 'id', [label '.id']);
entry.label = require_text_field(entry, 'label', [label '.label']);
entry.entry_type = optional_text_field(entry, 'entry_type', '');
entry.enabled = optional_logical_field(entry, 'enabled', true);
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
    entry.aggregate_ids = parse_aggregate_ids(entry.aggregate_ids, ...
        [label '.aggregate_ids']);
end
if ~isfield(entry, 'primary_particle_aggregate_ids') || ...
        isempty(entry.primary_particle_aggregate_ids)
    % Legacy profiles used aggregate_ids for both morphology and primary-
    % particle processing. Preserve that behavior when the separate field is
    % absent.
    entry.primary_particle_aggregate_ids = entry.aggregate_ids;
else
    entry.primary_particle_aggregate_ids = parse_aggregate_ids( ...
        entry.primary_particle_aggregate_ids, ...
        [label '.primary_particle_aggregate_ids']);
end
if ~isempty(entry.aggregate_ids) && ...
        any(~ismember(entry.primary_particle_aggregate_ids, ...
        entry.aggregate_ids))
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:PrimaryIdsOutsideMorphology', ...
        ['%s.primary_particle_aggregate_ids must be a subset of ', ...
        '%s.aggregate_ids.'], label, label);
end

% Canonical configs use style. Legacy plot.color/marker is upgraded here.
if isfield(entry, 'style') && ~isempty(entry.style)
    entry.style = validate_entry_style(entry.style, [label '.style']);
else
    legacy_plot = optional_struct_field(entry, 'plot');
    legacy_color = optional_text_field(legacy_plot, 'color', '');
    legacy_marker = optional_text_field(legacy_plot, 'marker', 'o');
    entry.style = legacy_style(entry.entry_type, legacy_color, legacy_marker);
end

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

function style = validate_entry_style(style, label)

if ~isstruct(style)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidStyle', ...
        '%s must be a struct.', label);
end
style.color = validate_hex_color(require_text_field(style, 'color', ...
    [label '.color']), [label '.color']);
style.secondary_color = validate_hex_color(require_text_field(style, ...
    'secondary_color', [label '.secondary_color']), ...
    [label '.secondary_color']);
style.box_face_color = optional_style_color(style, 'box_face_color', ...
    style.color, label);
style.median_color = optional_style_color(style, 'median_color', ...
    style.secondary_color, label);
style.marker = require_text_field(style, 'marker', [label '.marker']);
style.marker_size = require_positive_scalar_field(style, 'marker_size', ...
    [label '.marker_size']);
style.line_width = require_positive_scalar_field(style, 'line_width', ...
    [label '.line_width']);
style.box_face_alpha = require_unit_interval_field(style, ...
    'box_face_alpha', [label '.box_face_alpha']);
style.stack_colors = validate_stack_colors(style, label);

end

function stack_colors = validate_stack_colors(style, label)

if ~isfield(style, 'stack_colors') || ~isstruct(style.stack_colors)
    % Legacy profiles receive monochromatic stacks based on the condition color.
    stack_colors = struct('hybridity', {{style.color}}, ...
        'collapse', {{style.color}});
    return
end
stack_colors.hybridity = validate_color_list(style.stack_colors, ...
    'hybridity', [label '.stack_colors.hybridity']);
stack_colors.collapse = validate_color_list(style.stack_colors, ...
    'collapse', [label '.stack_colors.collapse']);

end

function color = optional_style_color(style, field_name, default_color, label)

color = default_color;
if isfield(style, field_name) && ~isempty(style.(field_name))
    color = validate_hex_color(style.(field_name), ...
        [label '.' field_name]);
end

end

function style = legacy_style(entry_type, explicit_color, marker)

% Backward compatibility only. Canonical configs never use this mapping.
switch lower(char(entry_type))
    case {'low_agglomeration', 'low_agglom', 'lal'}
        color = '#0072B2'; secondary = '#005A9C';
    case {'high_agglomeration', 'high_agglom', 'ex_agglomeration', 'exaglom'}
        color = '#D55E00'; secondary = '#C45100';
    case {'moderate_collapse', 'moderate_agglomeration', 'mod_collapse', 'hal'}
        color = '#009E73'; secondary = '#007A59';
    case {'extra_collapse', 'extreme_collapse', 'ex_collapse', 'excolaps'}
        color = '#CC79A7'; secondary = '#9E4F7F';
    otherwise
        color = '#0072B2'; secondary = '#005A9C';
end
if ~isempty(explicit_color)
    color = validate_hex_color(explicit_color, 'entries.plot.color');
    secondary = color;
end
style = struct('color', color, 'secondary_color', secondary, ...
    'box_face_color', color, 'median_color', secondary, ...
    'marker', marker, 'marker_size', 36, 'line_width', 1.5, ...
    'box_face_alpha', 0.25, 'stack_colors', struct( ...
    'hybridity', {{color}}, 'collapse', {{color}}));

end

function entries = validate_stack_palette_lengths(entries, frequency_bins)

% Each condition supplies one shade per category so every stack can remain
% within its condition hue while preserving the configured category order.
metric_names = {'hybridity','collapse'};
for i = 1:numel(entries)
    for j = 1:numel(metric_names)
        metric_name = metric_names{j};
        colors = entries(i).style.stack_colors.(metric_name);
        required_count = numel(frequency_bins.(metric_name));
        if isscalar(colors) && required_count > 1
            colors = repmat(colors, required_count, 1);
            entries(i).style.stack_colors.(metric_name) = colors;
        elseif numel(colors) ~= required_count
            error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidStackPalette', ...
                ['entries(%d).style.stack_colors.%s must contain one ', ...
                'color per configured bin (%d colors).'], ...
                i, metric_name, required_count);
        end
    end
end

end

function plots = validate_plots(cfg, entry_ids)

if ~isfield(cfg, 'plots') || ~isstruct(cfg.plots)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:MissingSection', ...
        'The TEM analysis config is missing the plots section.');
end
plots = cfg.plots;
plots.enabled = require_logical_field(plots, 'enabled', 'plots.enabled');
plots.visible = optional_text_field(plots, 'visible', 'on');
plots.export = require_logical_field(plots, 'export', 'plots.export');

% Validate the dual-format output schema while retaining legacy format and
% resolution fields for older custom profiles.
legacy_format = optional_text_field(plots, 'format', 'png');
plots.pdf = optional_logical_field(plots, 'pdf', strcmpi(legacy_format, 'pdf'));
plots.png = optional_logical_field(plots, 'png', ~strcmpi(legacy_format, 'pdf'));
legacy_resolution = optional_positive_integer_field(plots, 'resolution', 300);
plots.png_resolution = optional_positive_integer_field(plots, ...
    'png_resolution', legacy_resolution);
plots.png_source = lower(optional_text_field(plots, 'png_source', 'figure'));
if ~ismember(plots.png_source, {'figure','pdf'})
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidPNGSource', ...
        'plots.png_source must be figure or pdf.');
end
plots.save_figures = optional_logical_field(plots, 'save_figures', false);
if plots.export && ~plots.pdf && ~plots.png && ~plots.save_figures
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:NoFigureFormat', ...
        'At least one of plots.pdf, plots.png, or plots.save_figures must be enabled.');
end
if plots.export && plots.png && strcmp(plots.png_source, 'pdf') && ~plots.pdf
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:MissingPDFForPNG', ...
        'plots.pdf must be enabled when plots.png_source is pdf.');
end
plots.axis_line_width = optional_positive_scalar_field(plots, ...
    'axis_line_width', 1.0, 'plots.axis_line_width');
plots.axes = validate_axes_style(optional_struct_field(plots, 'axes'));
plots.boxplot = validate_boxplot_style(optional_struct_field(plots, ...
    'boxplot'));

if isfield(plots, 'font')
    plots.font = validate_font(plots.font);
else
    plots.font = struct('family', 'Segoe UI Semilight', 'fallback', 'Segoe UI', ...
        'axis_size', 13, 'label_size', 18, 'legend_size', 15, ...
        'title_size', 18, 'weight', 'normal', 'interpreter', 'tex');
end

if isfield(plots, 'reference')
    plots.reference = validate_reference(plots.reference);
else
    plots.reference = struct('color', '#7E2F8E', 'line_style', '-.', ...
        'line_width', 2.5);
end

if ~isfield(plots, 'frequency_bins') || ~isstruct(plots.frequency_bins)
    plots.frequency_bins = legacy_frequency_bins();
end
plots.frequency_bins.hybridity = validate_bins( ...
    plots.frequency_bins.hybridity, 'plots.frequency_bins.hybridity');
plots.frequency_bins.collapse = validate_bins( ...
    plots.frequency_bins.collapse, 'plots.frequency_bins.collapse');
validate_bin_coverage(plots.frequency_bins.hybridity, 1, Inf, ...
    'plots.frequency_bins.hybridity', true);
validate_bin_coverage(plots.frequency_bins.collapse, 0, 1, ...
    'plots.frequency_bins.collapse', false);

if ~isfield(plots, 'figures') || ~isstruct(plots.figures)
    plots.figures = legacy_figures(cfg.entries);
end
if ~isfield(plots.figures, 'subaggregate_count_distribution')
    plots.figures.subaggregate_count_distribution = ...
        legacy_subaggregate_count_distribution(cfg.entries);
end
required_figures = {'dpp_vs_da_manual_tem', ...
    'primary_particle_distributions', ...
    'appendix_a_primary_particle_distributions', ...
    'aggregate_metric_distributions', ...
    'subaggregate_count_distribution', ...
    'subaggregate_count_frequencies', ...
    'collapsed_subaggregate_frequencies', ...
    'dpp_vs_da_by_hybridity'};
for i = 1:numel(required_figures)
    id = required_figures{i};
    if ~isfield(plots.figures, id)
        error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:MissingFigureConfig', ...
            'plots.figures.%s is required.', id);
    end
    plots.figures.(id) = validate_figure_config( ...
        plots.figures.(id), id, entry_ids);
end
arrow_cfg = plots.figures.subaggregate_count_frequencies.arrows;
if arrow_cfg.target_bin_index > numel(plots.frequency_bins.hybridity)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidArrowBin', ...
        ['plots.figures.subaggregate_count_frequencies.' ...
        'arrows.target_bin_index exceeds the configured bin count.']);
end

end

function axes_style = validate_axes_style(axes_style)

axes_style.box = optional_logical_field(axes_style, 'box', true);
axes_style.layer = optional_text_field(axes_style, 'layer', 'top');
axes_style.tick_length = optional_numeric_vector(axes_style, 'tick_length', ...
    [0.02 0.02], 2, 'plots.axes.tick_length');
if any(axes_style.tick_length < 0)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidAxisStyle', ...
        'plots.axes.tick_length must contain nonnegative values.');
end
axes_style.tick_label_rotation = optional_nonnegative_scalar_field( ...
    axes_style, 'tick_label_rotation', 0, ...
    'plots.axes.tick_label_rotation');
if axes_style.tick_label_rotation > 360
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidAxisStyle', ...
        'plots.axes.tick_label_rotation must not exceed 360 degrees.');
end
axes_style.label_position_mode = optional_text_field(axes_style, ...
    'label_position_mode', 'auto');
if ~ismember(lower(axes_style.label_position_mode), {'auto','manual'})
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidAxisStyle', ...
        'plots.axes.label_position_mode must be auto or manual.');
end
axes_style.x_label_offset = optional_positive_scalar_field(axes_style, ...
    'x_label_offset', 0.055, 'plots.axes.x_label_offset');
axes_style.y_label_offset = optional_positive_scalar_field(axes_style, ...
    'y_label_offset', 0.075, 'plots.axes.y_label_offset');

end

function boxplot_style = validate_boxplot_style(boxplot_style)

boxplot_style.symbol = optional_text_field(boxplot_style, 'symbol', 'o');
boxplot_style.notch = optional_logical_field(boxplot_style, 'notch', true);
boxplot_style.median_min_line_width = optional_positive_scalar_field( ...
    boxplot_style, 'median_min_line_width', 2.0, ...
    'plots.boxplot.median_min_line_width');
boxplot_style.outlier_marker_size = optional_positive_scalar_field( ...
    boxplot_style, 'outlier_marker_size', 3, ...
    'plots.boxplot.outlier_marker_size');
boxplot_style.outlier_filled = optional_logical_field(boxplot_style, ...
    'outlier_filled', false);
boxplot_style.whisker_line_style = optional_text_field(boxplot_style, ...
    'whisker_line_style', '-');

end

function font = validate_font(font)

font.family = require_text_field(font, 'family', 'plots.font.family');
font.fallback = require_text_field(font, 'fallback', 'plots.font.fallback');
font.axis_size = require_positive_scalar_field(font, 'axis_size', ...
    'plots.font.axis_size');
font.label_size = require_positive_scalar_field(font, 'label_size', ...
    'plots.font.label_size');
font.legend_size = require_positive_scalar_field(font, 'legend_size', ...
    'plots.font.legend_size');
font.title_size = optional_positive_scalar_field(font, 'title_size', ...
    font.label_size, 'plots.font.title_size');
font.weight = require_text_field(font, 'weight', 'plots.font.weight');
if ~ismember(lower(font.weight), {'normal','bold'})
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidFontWeight', ...
        'plots.font.weight must be normal or bold.');
end
font.interpreter = require_text_field(font, 'interpreter', ...
    'plots.font.interpreter');
if ~ismember(lower(font.interpreter), {'tex','latex','none'})
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidInterpreter', ...
        'plots.font.interpreter must be tex, latex, or none.');
end

end

function reference = validate_reference(reference)

reference.color = validate_hex_color(require_text_field(reference, ...
    'color', 'plots.reference.color'), 'plots.reference.color');
reference.line_style = require_text_field(reference, 'line_style', ...
    'plots.reference.line_style');
reference.line_width = require_positive_scalar_field(reference, ...
    'line_width', 'plots.reference.line_width');

end

function colors = validate_color_list(src, field_name, label)

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:MissingColors', ...
        '%s must contain at least one color.', label);
end
colors = cellstr(src.(field_name));
for i = 1:numel(colors)
    colors{i} = validate_hex_color(colors{i}, sprintf('%s(%d)', label, i));
end

end

function bins = validate_bins(bins, label)

% Validate every interval independently, then compare all interval pairs so
% boundary inclusivity cannot create an overlap that depends on bin order.

if isempty(bins) || ~isstruct(bins)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidBins', ...
        '%s must be a non-empty struct array.', label);
end
ids = cell(numel(bins), 1);
for i = 1:numel(bins)
    item_label = sprintf('%s(%d)', label, i);
    bins(i).id = require_text_field(bins(i), 'id', [item_label '.id']);
    bins(i).label = require_text_field(bins(i), 'label', [item_label '.label']);
    bins(i).lower = optional_bound(bins(i), 'lower', -Inf, item_label);
    bins(i).upper = optional_bound(bins(i), 'upper', Inf, item_label);
    bins(i).include_lower = require_logical_field(bins(i), ...
        'include_lower', [item_label '.include_lower']);
    bins(i).include_upper = require_logical_field(bins(i), ...
        'include_upper', [item_label '.include_upper']);
    if bins(i).upper < bins(i).lower
        error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidBins', ...
            '%s upper bound must be >= lower bound.', item_label);
    end
    ids{i} = bins(i).id;
end
if numel(unique(ids)) ~= numel(ids)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:DuplicateBinId', ...
        '%s bin ids must be unique.', label);
end
for i = 1:numel(bins)
    for j = i + 1:numel(bins)
        if intervals_overlap(bins(i), bins(j))
            error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:OverlappingBins', ...
                '%s bins "%s" and "%s" overlap.', ...
                label, bins(i).id, bins(j).id);
        end
    end
end

end

function tf = intervals_overlap(a, b)

left = max(a.lower, b.lower);
right = min(a.upper, b.upper);
if left < right
    tf = true;
elseif left > right
    tf = false;
else
    in_a = (left > a.lower || a.include_lower) && ...
        (left < a.upper || a.include_upper);
    in_b = (left > b.lower || b.include_lower) && ...
        (left < b.upper || b.include_upper);
    tf = in_a && in_b;
end

end

function validate_bin_coverage(bins, domain_lower, domain_upper, label, ...
    discrete_domain)

[~, order] = sort([bins.lower]);
bins = bins(order);
if bins(1).lower > domain_lower || ...
        (bins(1).lower == domain_lower && ~bins(1).include_lower)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:IncompleteBins', ...
        '%s does not include its lower domain boundary.', label);
end
for i = 1:numel(bins) - 1
    if discrete_domain
        current_last = floor(bins(i).upper - double(~bins(i).include_upper));
        next_first = ceil(bins(i + 1).lower + ...
            double(~bins(i + 1).include_lower));
        gap_exists = next_first > current_last + 1;
    else
        gap_exists = bins(i).upper < bins(i + 1).lower || ...
            (bins(i).upper == bins(i + 1).lower && ...
            ~bins(i).include_upper && ~bins(i + 1).include_lower);
    end
    if gap_exists
        error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:IncompleteBins', ...
            '%s contains a gap between "%s" and "%s".', ...
            label, bins(i).id, bins(i + 1).id);
    end
end
if bins(end).upper < domain_upper || ...
        (isfinite(domain_upper) && bins(end).upper == domain_upper && ...
        ~bins(end).include_upper)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:IncompleteBins', ...
        '%s does not include its upper domain boundary.', label);
end

end

function value = optional_bound(src, field_name, default_value, label)

value = default_value;
if isfield(src, field_name) && ~isempty(src.(field_name))
    candidate = double(src.(field_name));
    % MATLAB decodes JSON null in numeric struct arrays as NaN.
    if isscalar(candidate) && isnan(candidate)
        return
    end
    value = candidate;
    if ~isscalar(value) || isnan(value)
        error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidBound', ...
            '%s.%s must be a numeric scalar or null.', label, field_name);
    end
end

end

function figure_cfg = validate_figure_config(figure_cfg, id, entry_ids)

% Every figure receives an independent ordered condition selection, complete
% axis definition, export name, layout style, and figure-specific options.

label = ['plots.figures.' id];
figure_cfg.enabled = require_logical_field(figure_cfg, 'enabled', ...
    [label '.enabled']);
figure_cfg.condition_ids = require_text_list(figure_cfg, ...
    'condition_ids', [label '.condition_ids']);
unknown = setdiff(figure_cfg.condition_ids, entry_ids);
if ~isempty(unknown)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:UnknownConditionId', ...
        '%s references unknown condition id(s): %s', ...
        label, strjoin(unknown, ', '));
end
figure_cfg.file_name = require_text_field(figure_cfg, 'file_name', ...
    [label '.file_name']);
figure_cfg.position = require_numeric_vector(figure_cfg, 'position', 4, ...
    [label '.position']);
if any(~isfinite(figure_cfg.position)) || any(figure_cfg.position(3:4) <= 0)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidPosition', ...
        '%s.position must have positive width and height.', label);
end
figure_cfg.x_scale = require_scale(figure_cfg, 'x_scale', label);
figure_cfg.y_scale = require_scale(figure_cfg, 'y_scale', label);
figure_cfg.x_limits = optional_axis_vector(figure_cfg, 'x_limits', 2, label);
figure_cfg.y_limits = optional_axis_vector(figure_cfg, 'y_limits', 2, label);
figure_cfg.x_ticks = optional_axis_vector(figure_cfg, 'x_ticks', [], label);
figure_cfg.y_ticks = optional_axis_vector(figure_cfg, 'y_ticks', [], label);
validate_axis_definition(figure_cfg.x_scale, figure_cfg.x_limits, ...
    figure_cfg.x_ticks, [label '.x']);
validate_axis_definition(figure_cfg.y_scale, figure_cfg.y_limits, ...
    figure_cfg.y_ticks, [label '.y']);
figure_cfg.x_minor_ticks = optional_axis_vector(figure_cfg, ...
    'x_minor_ticks', [], label);
validate_axis_definition(figure_cfg.x_scale, figure_cfg.x_limits, ...
    figure_cfg.x_minor_ticks, [label '.x_minor']);
figure_cfg.bar_width = optional_positive_scalar_field(figure_cfg, ...
    'bar_width', 0.4, [label '.bar_width']);
figure_cfg.condition_spacing = optional_positive_scalar_field(figure_cfg, ...
    'condition_spacing', 1.0, [label '.condition_spacing']);
figure_cfg.bar_edge_color = validate_named_or_hex_color( ...
    optional_text_field(figure_cfg, 'bar_edge_color', 'none'), ...
    [label '.bar_edge_color'], {'none'});
figure_cfg.bar_line_width = optional_nonnegative_scalar_field(figure_cfg, ...
    'bar_line_width', 0.5, [label '.bar_line_width']);
figure_cfg.layout_padding = optional_text_field(figure_cfg, ...
    'layout_padding', 'compact');
figure_cfg.tile_spacing = optional_text_field(figure_cfg, ...
    'tile_spacing', 'compact');
figure_cfg.box_width = optional_positive_scalar_field(figure_cfg, ...
    'box_width', 0.25, [label '.box_width']);
figure_cfg.ensemble_box_width = optional_positive_scalar_field(figure_cfg, ...
    'ensemble_box_width', 0.25, [label '.ensemble_box_width']);
figure_cfg.aggregate_box_width = optional_positive_scalar_field(figure_cfg, ...
    'aggregate_box_width', 0.3, [label '.aggregate_box_width']);
figure_cfg.kde_enabled = optional_logical_field(figure_cfg, ...
    'kde_enabled', false);
figure_cfg.kde_offset = optional_positive_scalar_field(figure_cfg, ...
    'kde_offset', 0.3, [label '.kde_offset']);
figure_cfg.kde_width = optional_positive_scalar_field(figure_cfg, ...
    'kde_width', 0.25, [label '.kde_width']);
figure_cfg.group_padding = optional_positive_scalar_field(figure_cfg, ...
    'group_padding', 0.3, [label '.group_padding']);
figure_cfg.legend_item_token_size = optional_numeric_vector(figure_cfg, ...
    'legend_item_token_size', [15 15], 2, ...
    [label '.legend_item_token_size']);
figure_cfg.axis_overrides = optional_struct_field(figure_cfg, 'axis_overrides');
figure_cfg.categories = optional_struct_array_field(figure_cfg, 'categories');

if isfield(figure_cfg, 'legend') && isstruct(figure_cfg.legend)
    legend_cfg = figure_cfg.legend;
else
    legend_cfg = struct();
end
figure_cfg.legend = struct( ...
    'location', optional_text_field(legend_cfg, 'location', 'best'), ...
    'orientation', optional_text_field(legend_cfg, 'orientation', 'vertical'), ...
    'columns', optional_positive_integer_field(legend_cfg, 'columns', 1), ...
    'box', optional_logical_field(legend_cfg, 'box', false), ...
    'interpreter', optional_text_field(legend_cfg, 'interpreter', 'none'));

if strcmp(id, 'subaggregate_count_distribution')
    if figure_cfg.enabled && numel(figure_cfg.condition_ids) ~= 2
        error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidMirroredConditions', ...
            ['%s.condition_ids must contain exactly two conditions when ', ...
            'the mirrored distribution is enabled.'], label);
    end
    if ~isempty(figure_cfg.x_limits) && ...
            (figure_cfg.x_limits(1) >= 0 || figure_cfg.x_limits(2) <= 0)
        error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidMirroredLimits', ...
            '%s.x_limits must span zero.', label);
    end
    figure_cfg.grid = validate_mirrored_grid( ...
        optional_struct_field(figure_cfg, 'grid'), label);
    figure_cfg.zero_line = validate_zero_line( ...
        optional_struct_field(figure_cfg, 'zero_line'), label);
    figure_cfg.stem_line_width = optional_positive_scalar_field(figure_cfg, ...
        'stem_line_width', 1.5, [label '.stem_line_width']);
    figure_cfg.marker_size = optional_positive_scalar_field(figure_cfg, ...
        'marker_size', 7, [label '.marker_size']);
    figure_cfg.marker_edge_width = optional_positive_scalar_field( ...
        figure_cfg, 'marker_edge_width', 1.5, ...
        [label '.marker_edge_width']);
    figure_cfg.marker_face_color = validate_named_or_hex_color( ...
        optional_text_field(figure_cfg, 'marker_face_color', 'none'), ...
        [label '.marker_face_color'], {'none'});
else
    figure_cfg.grid = struct();
    figure_cfg.zero_line = struct();
end

if strcmp(id, 'dpp_vs_da_by_hybridity')
    if isempty(figure_cfg.categories)
        error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:MissingCategories', ...
            '%s.categories must be configured.', label);
    end
    if ~isfield(figure_cfg.categories, 'line_width')
        [figure_cfg.categories.line_width] = deal(1.5);
    end
    for i = 1:numel(figure_cfg.categories)
        category_label = sprintf('%s.categories(%d)', label, i);
        c = figure_cfg.categories(i);
        c.id = require_text_field(c, 'id', [category_label '.id']);
        c.label = require_text_field(c, 'label', [category_label '.label']);
        c.lower = optional_bound(c, 'lower', -Inf, category_label);
        c.upper = optional_bound(c, 'upper', Inf, category_label);
        c.include_lower = require_logical_field(c, 'include_lower', ...
            [category_label '.include_lower']);
        c.include_upper = require_logical_field(c, 'include_upper', ...
            [category_label '.include_upper']);
        c.color = validate_hex_color(require_text_field(c, 'color', ...
            [category_label '.color']), [category_label '.color']);
        c.marker = require_text_field(c, 'marker', [category_label '.marker']);
        c.marker_size = require_positive_scalar_field(c, 'marker_size', ...
            [category_label '.marker_size']);
        c.line_width = require_positive_scalar_field(c, 'line_width', ...
            [category_label '.line_width']);
        figure_cfg.categories(i) = c;
    end
    validate_bins(figure_cfg.categories, [label '.categories']);
end

if strcmp(id, 'subaggregate_count_frequencies')
    figure_cfg.arrows = validate_arrow_config(figure_cfg, label, entry_ids);
else
    figure_cfg.arrows = struct();
end

end

function grid_cfg = validate_mirrored_grid(grid_cfg, label)

grid_cfg.x_major = optional_logical_field(grid_cfg, 'x_major', true);
grid_cfg.x_minor = optional_logical_field(grid_cfg, 'x_minor', true);
grid_cfg.y_major = optional_logical_field(grid_cfg, 'y_major', false);
grid_cfg.y_minor = optional_logical_field(grid_cfg, 'y_minor', false);
grid_cfg.major_color = validate_hex_color(optional_text_field(grid_cfg, ...
    'major_color', '#D9D9D9'), [label '.grid.major_color']);
grid_cfg.minor_color = validate_hex_color(optional_text_field(grid_cfg, ...
    'minor_color', '#EEEEEE'), [label '.grid.minor_color']);
grid_cfg.major_alpha = optional_unit_interval_field(grid_cfg, ...
    'major_alpha', 0.65, [label '.grid.major_alpha']);
grid_cfg.minor_alpha = optional_unit_interval_field(grid_cfg, ...
    'minor_alpha', 0.50, [label '.grid.minor_alpha']);
grid_cfg.major_line_style = optional_text_field(grid_cfg, ...
    'major_line_style', '-');
grid_cfg.minor_line_style = optional_text_field(grid_cfg, ...
    'minor_line_style', ':');

end

function zero_line = validate_zero_line(zero_line, label)

zero_line.color = validate_hex_color(optional_text_field(zero_line, ...
    'color', '#4D4D4D'), [label '.zero_line.color']);
zero_line.line_width = optional_positive_scalar_field(zero_line, ...
    'line_width', 1.0, [label '.zero_line.line_width']);
zero_line.line_style = optional_text_field(zero_line, ...
    'line_style', '-');

end

function arrow_cfg = validate_arrow_config(figure_cfg, label, entry_ids)

if isfield(figure_cfg, 'arrows') && isstruct(figure_cfg.arrows)
    arrow_cfg = figure_cfg.arrows;
else
    arrow_cfg = struct();
end
arrow_cfg.enabled = optional_logical_field(arrow_cfg, 'enabled', true);
if isfield(arrow_cfg, 'condition_ids') && ~isempty(arrow_cfg.condition_ids)
    arrow_cfg.condition_ids = require_text_list(arrow_cfg, ...
        'condition_ids', [label '.arrows.condition_ids']);
else
    count = min(2, numel(figure_cfg.condition_ids));
    arrow_cfg.condition_ids = figure_cfg.condition_ids(1:count);
end
unknown = setdiff(arrow_cfg.condition_ids, entry_ids);
if ~isempty(unknown)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:UnknownConditionId', ...
        '%s.arrows references unknown condition id(s): %s', ...
        label, strjoin(unknown, ', '));
end
outside_figure = setdiff(arrow_cfg.condition_ids, figure_cfg.condition_ids);
if ~isempty(outside_figure)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidArrowCondition', ...
        '%s.arrows condition ids must be selected by the figure.', ...
        label);
end
arrow_cfg.target_bin_index = optional_positive_integer_field(arrow_cfg, ...
    'target_bin_index', 1);
arrow_cfg.tail_fraction = optional_unit_interval_field(arrow_cfg, ...
    'tail_fraction', 1.0, [label '.arrows.tail_fraction']);
arrow_cfg.shaft_length = optional_positive_scalar_field(arrow_cfg, ...
    'shaft_length', 10, [label '.arrows.shaft_length']);
arrow_cfg.minimum_head_clearance = optional_nonnegative_scalar_field( ...
    arrow_cfg, 'minimum_head_clearance', 2, ...
    [label '.arrows.minimum_head_clearance']);
arrow_cfg.maximum_tail_y = optional_positive_scalar_field(arrow_cfg, ...
    'maximum_tail_y', 98, [label '.arrows.maximum_tail_y']);
arrow_cfg.line_width = optional_positive_scalar_field(arrow_cfg, ...
    'line_width', 2.0, [label '.arrows.line_width']);
arrow_cfg.color = validate_hex_color(optional_text_field(arrow_cfg, ...
    'color', '#000000'), [label '.arrows.color']);
arrow_cfg.head_marker = optional_text_field(arrow_cfg, ...
    'head_marker', 'v');
arrow_cfg.head_marker_size = optional_positive_scalar_field(arrow_cfg, ...
    'head_marker_size', 8, [label '.arrows.head_marker_size']);

end

function validate_axis_definition(scale, limits, ticks, label)

if ~isempty(ticks) && any(diff(ticks) <= 0)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidAxisVector', ...
        '%s_ticks must be strictly increasing.', label);
end
if strcmpi(scale, 'log') && (any(limits <= 0) || any(ticks <= 0))
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidAxisVector', ...
        '%s limits and ticks must be positive on a log scale.', label);
end
if ~isempty(limits) && ~isempty(ticks) && ...
        (ticks(1) < limits(1) || ticks(end) > limits(2))
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidAxisVector', ...
        '%s ticks must remain within the configured limits.', label);
end

end

function value = require_scale(src, field_name, label)

value = optional_text_field(src, field_name, 'linear');
if ~ismember(lower(value), {'linear','log'})
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidScale', ...
        '%s.%s must be linear or log.', label, field_name);
end

end

function values = optional_axis_vector(src, field_name, expected_length, label)

values = [];
if ~isfield(src, field_name) || isempty(src.(field_name))
    return
end
values = double(src.(field_name)(:)).';
if any(~isfinite(values)) || (~isempty(expected_length) && ...
        numel(values) ~= expected_length)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidAxisVector', ...
        '%s.%s has invalid values or length.', label, field_name);
end
if numel(values) == 2 && values(2) <= values(1)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidAxisVector', ...
        '%s.%s must be increasing.', label, field_name);
end

end

function values = optional_numeric_vector(src, field_name, default_values, ...
    expected_length, label)

values = default_values;
if isfield(src, field_name) && ~isempty(src.(field_name))
    values = double(src.(field_name)(:)).';
    if numel(values) ~= expected_length || any(~isfinite(values))
        error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidVector', ...
            '%s must contain %d finite numeric values.', ...
            label, expected_length);
    end
end

end

function out = require_numeric_vector(src, field_name, length_value, label)

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:MissingVector', ...
        '%s is required.', label);
end
out = double(src.(field_name)(:)).';
if numel(out) ~= length_value
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidVector', ...
        '%s must contain %d values.', label, length_value);
end

end

function values = require_text_list(src, field_name, label)

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:MissingTextList', ...
        '%s must be a non-empty text array.', label);
end
values = cellstr(src.(field_name));
if numel(unique(values)) ~= numel(values)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:DuplicateConditionId', ...
        '%s must not contain duplicates.', label);
end

end

function bins = legacy_frequency_bins()

bins.hybridity = struct( ...
    'id', {'eq_1','eq_2','three_to_five','six_to_ten','gt_10'}, ...
    'label', {'{\it N}_{\rm subagg} = 1','{\it N}_{\rm subagg} = 2', ...
        '3 \leq {\it N}_{\rm subagg} \leq 5', ...
        '6 \leq {\it N}_{\rm subagg} \leq 10', ...
        '{\it N}_{\rm subagg} > 10'}, ...
    'lower', {1,2,3,6,10}, 'upper', {1,2,5,10,[]}, ...
    'include_lower', {true,true,true,true,false}, ...
    'include_upper', {true,true,true,true,false});
bins.collapse = struct( ...
    'id', {'eq_0','zero_to_low','low_to_high','high_to_one','eq_1'}, ...
    'label', {'{\it f}_{\rm col,subagg} = 0', ...
        '0 < {\it f}_{\rm col,subagg} < 0.33', ...
        '0.33 \leq {\it f}_{\rm col,subagg} \leq 0.67', ...
        '0.67 < {\it f}_{\rm col,subagg} < 1', ...
        '{\it f}_{\rm col,subagg} = 1'}, ...
    'lower', {0,0,0.33,0.67,1}, 'upper', {0,0.33,0.67,1,1}, ...
    'include_lower', {true,false,true,false,true}, ...
    'include_upper', {true,false,true,false,true});

end

function figures = legacy_figures(entries)

% Construct a complete in-memory figure schema only for legacy profiles that
% predate per-figure configuration. Canonical profiles define these fields.

ids = {entries([entries.enabled]).id};
if isempty(ids)
    ids = {entries.id};
end
base = struct('enabled', true, 'condition_ids', {ids}, ...
    'file_name', '', 'position', [100 100 900 600], ...
    'x_scale', 'linear', 'y_scale', 'linear', ...
    'x_limits', [], 'y_limits', [], 'x_ticks', [], 'y_ticks', [], ...
    'legend', struct('location','best','orientation','vertical', ...
        'columns',1,'box',false,'interpreter','none'), ...
    'axis_overrides', struct());
names = {'dpp_vs_da_manual_tem','primary_particle_distributions', ...
    'appendix_a_primary_particle_distributions', ...
    'aggregate_metric_distributions','subaggregate_count_frequencies', ...
    'collapsed_subaggregate_frequencies','dpp_vs_da_by_hybridity'};
positions = {[50 50 900 500],[100 0 850 1800],[100 0 850 1800], ...
    [150 50 1050 900],[200 150 900 500],[200 150 900 500], ...
    [250 200 900 500]};
for i = 1:numel(names)
    item = base;
    item.file_name = names{i};
    item.position = positions{i};
    figures.(names{i}) = item;
end
figures.subaggregate_count_distribution = ...
    legacy_subaggregate_count_distribution(entries);
figures.appendix_a_primary_particle_distributions.enabled = false;
figures.subaggregate_count_frequencies.legend.interpreter = 'tex';
figures.collapsed_subaggregate_frequencies.legend.interpreter = 'tex';
figures.dpp_vs_da_by_hybridity.legend.interpreter = 'tex';
figures.dpp_vs_da_by_hybridity.categories = struct( ...
    'id', {'eq_1','eq_2','three_to_five','gt_5'}, ...
    'label', {'{\it N}_{\rm subagg} = 1', ...
        '{\it N}_{\rm subagg} = 2', ...
        '3 \leq {\it N}_{\rm subagg} \leq 5', ...
        '{\it N}_{\rm subagg} > 5'}, ...
    'lower', {1,2,3,5}, 'upper', {1,2,5,[]}, ...
    'include_lower', {true,true,true,false}, ...
    'include_upper', {true,true,true,false}, ...
    'color', {'#0072B2','#E69F00','#009E73','#CC79A7'}, ...
    'marker', {'^','s','h','o'}, 'marker_size', {36,36,48,36}, ...
    'line_width', {1.5,1.5,1.5,1.5});

end

function figure_cfg = legacy_subaggregate_count_distribution(entries)

% Older custom profiles remain valid while the mirrored figure stays opt-in
% until its two comparison conditions are selected explicitly.
ids = {entries([entries.enabled]).id};
if isempty(ids)
    ids = {entries.id};
end
ids = ids(1:min(2, numel(ids)));
figure_cfg = struct( ...
    'enabled', false, ...
    'condition_ids', {ids}, ...
    'file_name', 'subaggregate_count_distribution', ...
    'position', [200 0 1000 1800], ...
    'x_scale', 'linear', ...
    'y_scale', 'log', ...
    'x_limits', [-100 100], ...
    'y_limits', [0.85 120], ...
    'x_ticks', -100:20:100, ...
    'x_minor_ticks', [], ...
    'y_ticks', [1 2 3 4 5 10:10:100], ...
    'stem_line_width', 1.5, ...
    'marker_size', 4, ...
    'marker_edge_width', 1.2, ...
    'marker_face_color', 'none', ...
    'legend_item_token_size', [18 12], ...
    'legend', struct('location', 'eastoutside', ...
        'orientation', 'vertical', 'columns', 1, 'box', false, ...
        'interpreter', 'tex'), ...
    'grid', struct('x_major', false, 'x_minor', false, ...
        'y_major', false, 'y_minor', false), ...
    'zero_line', struct('color', '#4D4D4D', 'line_width', 1.0, ...
        'line_style', '-'));

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
            [label '.start']);
        end_id = require_positive_integer_field(src, 'end', ...
            [label '.end']);
        step_id = optional_positive_integer_field(src, 'step', 1);
        if end_id < start_id
            error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidAggregateIdRange', ...
                '%s.end must be greater than or equal to start.', ...
                label);
        end
        ids = (start_id:step_id:end_id).';
    end
else
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidAggregateIds', ...
        '%s must be a numeric array or range struct.', label);
end
if isempty(ids) || any(~isfinite(ids)) || any(ids <= 0) || any(round(ids) ~= ids)
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidAggregateIds', ...
        '%s must contain positive integer ids.', label);
end
ids = unique(ids(:), 'stable');

end

function out = require_struct_field(src, field_name, label)

% Required-struct helpers centralize type checks so loader errors identify the
% exact configuration field instead of failing later in the analysis script.

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

function out = optional_struct_array_field(src, field_name)

out = struct([]);
if isfield(src, field_name) && ~isempty(src.(field_name))
    if ~isstruct(src.(field_name))
        error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidStruct', ...
            'The TEM analysis config field %s must be a struct array.', field_name);
    end
    out = src.(field_name);
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

function values = require_scalar_field(src, field_name, label)

if ~isfield(src, field_name) || isempty(src.(field_name))
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:MissingScalar', ...
        'The TEM analysis config is missing %s.', label);
end
values = double(src.(field_name));
if ~isscalar(values) || ~isfinite(values)
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

function out = optional_positive_scalar_field(src, field_name, default_value, label)

out = default_value;
if isfield(src, field_name) && ~isempty(src.(field_name))
    out = require_positive_scalar_field(src, field_name, label);
end

end

function out = optional_nonnegative_scalar_field(src, field_name, ...
    default_value, label)

out = default_value;
if isfield(src, field_name) && ~isempty(src.(field_name))
    out = require_scalar_field(src, field_name, label);
    if out < 0
        error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidNonnegativeScalar', ...
            'The TEM analysis config field %s must be nonnegative.', label);
    end
end

end

function out = optional_unit_interval_field(src, field_name, ...
    default_value, label)

out = default_value;
if isfield(src, field_name) && ~isempty(src.(field_name))
    out = require_unit_interval_field(src, field_name, label);
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

function out = require_unit_interval_field(src, field_name, label)

out = require_scalar_field(src, field_name, label);
if out < 0 || out > 1
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidUnitInterval', ...
        '%s must be between 0 and 1.', label);
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

function value = validate_hex_color(value, label)

value = upper(char(value));
if isempty(regexp(value, '^#[0-9A-F]{6}$', 'once'))
    error('PFAL:LOAD_MAIN_TEM_ANALYSIS_CONFIG:InvalidColor', ...
        '%s must be a six-digit hexadecimal color.', label);
end

end

function value = validate_named_or_hex_color(value, label, allowed_names)

% Permit explicitly supported MATLAB color controls alongside hexadecimal
% presentation colors without accepting ambiguous color abbreviations.
value = char(value);
if any(strcmpi(value, allowed_names))
    value = lower(value);
    return
end
value = validate_hex_color(value, label);

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
