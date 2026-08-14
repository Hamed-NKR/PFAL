clc
clear
close all
warning('off')

%% Load config

cfg_tem = UTILS.LOAD_MAIN_TEM_ANALYSIS_CONFIG;

if ~isempty(cfg_tem.analysis.random_seed)
    rng(cfg_tem.analysis.random_seed)
end

if ~isfolder(cfg_tem.outputs.data_root)
    mkdir(cfg_tem.outputs.data_root)
end
if cfg_tem.plots.enabled && cfg_tem.plots.export && ...
        ~isfolder(cfg_tem.outputs.results_root)
    mkdir(cfg_tem.outputs.results_root)
end

%% Process enabled TEM entries

% Disabled entries remain in the saved configuration for future condition
% expansion, but they do not participate in data loading or analysis.
active_entries = cfg_tem.entries([cfg_tem.entries.enabled]);
n_entries = numel(active_entries);
entry_result_cells = cell(n_entries, 1);
aggregate_tables = cell(n_entries, 1);
morphology_tables = cell(n_entries, 1);
primary_tables = cell(n_entries, 1);
summary_tables = cell(n_entries, 1);
model_input_tables = cell(n_entries, 1);

for i = 1:n_entries
    entry_result_cells{i} = process_tem_entry(active_entries(i), cfg_tem.analysis);
    aggregate_tables{i} = entry_result_cells{i}.aggregate_table;
    morphology_tables{i} = entry_result_cells{i}.morphology_table;
    primary_tables{i} = entry_result_cells{i}.primary_table;
    summary_tables{i} = entry_result_cells{i}.summary_table;
    model_input_tables{i} = entry_result_cells{i}.model_inputs;
end

entry_results = vertcat(entry_result_cells{:});
aggregate_table = vertcat_nonempty(aggregate_tables);
morphology_table = vertcat_nonempty(morphology_tables);
primary_table = vertcat_nonempty(primary_tables);
summary_table = vertcat_nonempty(summary_tables);
model_inputs = vertcat_nonempty(model_input_tables);

% Build figure-facing summaries separately from the scalar model inputs so
% the plotted quantities remain available for auditing and reuse.
plot_artifacts = build_plot_artifacts(entry_results, cfg_tem);
ensemble_primary_summary = plot_artifacts.ensemble_primary_summary;
per_aggregate_dpp_summary = plot_artifacts.per_aggregate_dpp_summary;
per_aggregate_sigmapp_summary = plot_artifacts.per_aggregate_sigmapp_summary;
hybridity_frequency_summary = plot_artifacts.hybridity_frequency_summary;
collapse_frequency_summary = plot_artifacts.collapse_frequency_summary;
subaggregate_count_distribution_summary = ...
    plot_artifacts.subaggregate_count_distribution_summary;
plot_diagnostics = plot_artifacts.plot_diagnostics;

%% Export machine-readable outputs

metadata = struct( ...
    'generated_at', char(datetime('now')), ...
    'config_file', cfg_tem.config_file, ...
    'entry_count', n_entries, ...
    'configured_entry_count', numel(cfg_tem.entries));

save(fullfile(cfg_tem.outputs.data_root, cfg_tem.outputs.mat_file), ...
    'cfg_tem', 'metadata', 'entry_results', 'aggregate_table', ...
    'morphology_table', 'primary_table', 'summary_table', 'model_inputs', ...
    'plot_artifacts', ...
    'ensemble_primary_summary', 'per_aggregate_dpp_summary', ...
    'per_aggregate_sigmapp_summary', 'hybridity_frequency_summary', ...
    'collapse_frequency_summary', ...
    'subaggregate_count_distribution_summary', 'plot_diagnostics', '-v7.3')

writetable(summary_table, fullfile(cfg_tem.outputs.data_root, ...
    cfg_tem.outputs.summary_csv))
writetable(model_inputs, fullfile(cfg_tem.outputs.data_root, ...
    cfg_tem.outputs.model_inputs_csv))
write_optional_table(ensemble_primary_summary, fullfile( ...
    cfg_tem.outputs.data_root, 'primary_ensemble_summary.csv'))
write_optional_table(per_aggregate_dpp_summary, fullfile( ...
    cfg_tem.outputs.data_root, 'per_aggregate_dpp_summary.csv'))
write_optional_table(per_aggregate_sigmapp_summary, fullfile( ...
    cfg_tem.outputs.data_root, 'per_aggregate_sigmapp_summary.csv'))
write_optional_table(hybridity_frequency_summary, fullfile( ...
    cfg_tem.outputs.data_root, 'hybridity_frequency_summary.csv'))
write_optional_table(collapse_frequency_summary, fullfile( ...
    cfg_tem.outputs.data_root, 'collapse_frequency_summary.csv'))
distribution_cfg = ...
    cfg_tem.plots.figures.subaggregate_count_distribution;
write_optional_table(subaggregate_count_distribution_summary, fullfile( ...
    cfg_tem.outputs.data_root, ...
    [distribution_cfg.file_name '_summary.csv']))
write_optional_table(plot_diagnostics, fullfile( ...
    cfg_tem.outputs.data_root, 'plot_diagnostics.csv'))

%% Plot outputs

if cfg_tem.plots.enabled
    plot_tem_analysis(entry_results, cfg_tem, plot_artifacts)
    if strcmpi(cfg_tem.plots.visible, 'off')
        close all force
    end
end

%% Local functions

function entry_result = process_tem_entry(entry, analysis_cfg)

if ~isfile(entry.aggregate.file)
    error('PFAL:main_tem_analysis_v1:MissingAggregateFile', ...
        'Aggregate MAT file not found for entry "%s": %s', ...
        entry.id, entry.aggregate.file);
end

loaded = load(entry.aggregate.file, entry.aggregate.variable);
if ~isfield(loaded, entry.aggregate.variable)
    error('PFAL:main_tem_analysis_v1:MissingAggregateVariable', ...
        'Variable "%s" was not found in %s.', ...
        entry.aggregate.variable, entry.aggregate.file);
end
Aggs = loaded.(entry.aggregate.variable);

if isempty(entry.aggregate_ids)
    morphology_ids = (1:numel(Aggs)).';
else
    morphology_ids = entry.aggregate_ids(:);
end
if isempty(entry.primary_particle_aggregate_ids)
    primary_particle_ids = (1:numel(Aggs)).';
else
    primary_particle_ids = entry.primary_particle_aggregate_ids(:);
end

if any(morphology_ids > numel(Aggs))
    error('PFAL:main_tem_analysis_v1:AggregateIdOutOfRange', ...
        ['Entry "%s" includes morphology aggregate ids larger than ', ...
        'numel(Aggs) = %d.'], ...
        entry.id, numel(Aggs));
end
if any(primary_particle_ids > numel(Aggs))
    error('PFAL:main_tem_analysis_v1:PrimaryParticleAggregateIdOutOfRange', ...
        ['Entry "%s" includes primary-particle aggregate ids larger than ', ...
        'numel(Aggs) = %d.'], ...
        entry.id, numel(Aggs));
end

aggregate_rows = table();
primary_rows = table();
morphology_rows = build_morphology_table(Aggs, morphology_ids, entry);

for j = 1:numel(primary_particle_ids)
    agg_id = primary_particle_ids(j);
    csv_file = fullfile(entry.primary_particles.folder, ...
        sprintf(entry.primary_particles.file_pattern, agg_id));

    if ~isfile(csv_file)
        warning('PFAL:main_tem_analysis_v1:MissingPrimaryCsv', ...
            'Primary-particle CSV not found for entry "%s", aggregate %d: %s', ...
            entry.id, agg_id, csv_file);
        continue
    end

    pp_table = readtable(csv_file);
    if ~ismember(entry.primary_particles.area_column, pp_table.Properties.VariableNames)
        error('PFAL:main_tem_analysis_v1:MissingAreaColumn', ...
            'CSV file %s is missing the "%s" area column.', ...
            csv_file, entry.primary_particles.area_column);
    end

    area_values = pp_table.(entry.primary_particles.area_column);
    first_row = analysis_cfg.primary_area_start_row;
    if numel(area_values) < first_row
        warning('PFAL:main_tem_analysis_v1:EmptyPrimaryCsv', ...
            'CSV file %s has no primary-particle rows after row %d.', ...
            csv_file, first_row - 1);
        continue
    end

    area_values = double(area_values(first_row:end));
    area_values = area_values(isfinite(area_values) & area_values > 0);
    if isempty(area_values)
        warning('PFAL:main_tem_analysis_v1:NoPositivePrimaryArea', ...
            'CSV file %s has no positive primary-particle areas.', csv_file);
        continue
    end

    % Convert ImageJ primary-particle projected areas to equivalent-area
    % diameters.  The first CSV row is intentionally skipped by config
    % because the manual ImageJ exports include an aggregate/background row
    % before the individual primary particles.
    dpp_values = sqrt(4 * area_values / pi);
    da_value = require_agg_numeric_field(Aggs, agg_id, 'da');

    % Coverage is the ratio of total manually sized primary-particle area to
    % aggregate projected area. The smooth logistic switch applies a stronger
    % inverse-area correction at low coverage and approaches unweighted
    % statistics continuously at high coverage.
    coverage = sum(dpp_values .^ 2) / (da_value ^ 2);
    alpha_value = 1 - 1 / (1 + exp( ...
        -(coverage - analysis_cfg.coverage_threshold) / ...
        analysis_cfg.logistic_bandwidth));

    % Per-aggregate weights estimate the geometric mean and GSD of primary
    % particles within this one aggregate.  Ensemble weights below include
    % parent aggregate area because the ensemble distribution represents a
    % condition-level population rather than equal-weight manual measurements.
    per_aggregate_weights = (1 ./ (dpp_values .^ 2)) .^ alpha_value;
    [dbarpp, sigmapp, ci_dbarpp] = UTILS.GEOMSTATS( ...
        dpp_values, per_aggregate_weights);

    npp = numel(dpp_values);
    parent_da = repmat(da_value, npp, 1);
    parent_cov = repmat(coverage, npp, 1);
    parent_alpha = repmat(alpha_value, npp, 1);
    ensemble_weight_soft = ((parent_da .^ 2) ./ ...
        (dpp_values .^ 2)) .^ parent_alpha;
    ensemble_weight_hard = ones(npp, 1);
    hard_mask = parent_cov < analysis_cfg.coverage_threshold;
    ensemble_weight_hard(hard_mask) = ...
        (parent_da(hard_mask) .^ 2) ./ (dpp_values(hard_mask) .^ 2);

    primary_rows = [primary_rows; table( ...
        repmat({entry.id}, npp, 1), ...
        repmat({entry.label}, npp, 1), ...
        repmat(agg_id, npp, 1), ...
        (1:npp).', ...
        area_values(:), ...
        dpp_values(:), ...
        parent_da(:), ...
        parent_cov(:), ...
        ensemble_weight_hard(:), ...
        ensemble_weight_soft(:), ...
        repmat({csv_file}, npp, 1), ...
        'VariableNames', {'entry_id', 'entry_label', 'aggregate_id', ...
        'primary_index', 'area', 'dpp_nm', 'parent_da_nm', ...
        'coverage', 'weight_hard', 'weight_soft', 'csv_file'})]; %#ok<AGROW>

    aggregate_rows = [aggregate_rows; table( ...
        {entry.id}, ...
        {entry.label}, ...
        agg_id, ...
        da_value, ...
        npp, ...
        dbarpp, ...
        sigmapp, ...
        ci_dbarpp(1), ...
        ci_dbarpp(2), ...
        coverage, ...
        alpha_value, ...
        optional_agg_numeric_field(Aggs, agg_id, 'n_subagg'), ...
        optional_agg_numeric_field(Aggs, agg_id, 'n_colaps'), ...
        optional_agg_numeric_field(Aggs, agg_id, 'ca'), ...
        optional_agg_numeric_field(Aggs, agg_id, 'zbar_opt'), ...
        optional_agg_numeric_field(Aggs, agg_id, 'sbar_opt'), ...
        {csv_file}, ...
        'VariableNames', {'entry_id', 'entry_label', 'aggregate_id', ...
        'da_nm', 'npp_manual', 'dbarpp_nm', 'sigmapp', ...
        'dbarpp_ci95_low_nm', 'dbarpp_ci95_high_nm', 'coverage', ...
        'coverage_weight_alpha', 'n_subagg', 'n_colaps', 'ca', ...
        'zbar_opt', 'sbar_opt', 'csv_file'})]; %#ok<AGROW>
end

if isempty(aggregate_rows)
    error('PFAL:main_tem_analysis_v1:NoAggregateRows', ...
        'Entry "%s" did not produce any aggregate rows.', entry.id);
end

summary_table = summarize_entry(entry, aggregate_rows, morphology_rows, ...
    primary_rows, analysis_cfg);
model_inputs = build_model_inputs(entry, aggregate_rows);

entry_result = struct();
entry_result.entry = entry;
entry_result.aggregate_table = aggregate_rows;
entry_result.morphology_table = morphology_rows;
entry_result.primary_table = primary_rows;
entry_result.summary_table = summary_table;
entry_result.model_inputs = model_inputs;

end

function morphology_table = build_morphology_table(Aggs, aggregate_ids, entry)

% Morphology records are loaded directly from the aggregate MAT file and do
% not depend on the availability of manually sized primary-particle CSVs.
n_rows = numel(aggregate_ids);
entry_id = repmat({entry.id}, n_rows, 1);
entry_label = repmat({entry.label}, n_rows, 1);
aggregate_id = aggregate_ids(:);
da_nm = NaN(n_rows, 1);
n_subagg = NaN(n_rows, 1);
n_colaps = NaN(n_rows, 1);
ca = NaN(n_rows, 1);
zbar_opt = NaN(n_rows, 1);
sbar_opt = NaN(n_rows, 1);
for i = 1:n_rows
    agg_id = aggregate_id(i);
    da_nm(i) = optional_agg_numeric_field(Aggs, agg_id, 'da');
    n_subagg(i) = optional_agg_numeric_field(Aggs, agg_id, 'n_subagg');
    n_colaps(i) = optional_agg_numeric_field(Aggs, agg_id, 'n_colaps');
    ca(i) = optional_agg_numeric_field(Aggs, agg_id, 'ca');
    zbar_opt(i) = optional_agg_numeric_field(Aggs, agg_id, 'zbar_opt');
    sbar_opt(i) = optional_agg_numeric_field(Aggs, agg_id, 'sbar_opt');
end
morphology_table = table(entry_id, entry_label, aggregate_id, da_nm, ...
    n_subagg, n_colaps, ca, zbar_opt, sbar_opt);

end

function summary_table = summarize_entry(entry, aggregate_table, ...
    morphology_table, primary_table, analysis_cfg)

[GM_dpp_ens_unweighted, GSD_dpp_ens_unweighted, CI_dpp_ens_unweighted] = ...
    UTILS.GEOMSTATS(primary_table.dpp_nm);
[GM_dpp_ens_hard, GSD_dpp_ens_hard, CI_dpp_ens_hard] = ...
    UTILS.GEOMSTATS(primary_table.dpp_nm, primary_table.weight_hard);
[GM_dpp_ens_soft, GSD_dpp_ens_soft, CI_dpp_ens_soft] = ...
    UTILS.GEOMSTATS(primary_table.dpp_nm, primary_table.weight_soft);

[GM_dbarpp, GSD_dbarpp, CI_dbarpp] = ...
    UTILS.GEOMSTATS(aggregate_table.dbarpp_nm);
[GM_sigmapp, GSD_sigmapp, CI_sigmapp] = ...
    UTILS.GEOMSTATS(aggregate_table.sigmapp);
[GM_da, GSD_da, CI_da] = UTILS.GEOMSTATS(aggregate_table.da_nm);

summary_table = table( ...
    {entry.id}, ...
    {entry.label}, ...
    height(aggregate_table), ...
    height(morphology_table), ...
    height(primary_table), ...
    analysis_cfg.coverage_threshold, ...
    analysis_cfg.logistic_bandwidth, ...
    GM_dpp_ens_unweighted, GSD_dpp_ens_unweighted, ...
    CI_dpp_ens_unweighted(1), CI_dpp_ens_unweighted(2), ...
    GM_dpp_ens_hard, GSD_dpp_ens_hard, ...
    CI_dpp_ens_hard(1), CI_dpp_ens_hard(2), ...
    GM_dpp_ens_soft, GSD_dpp_ens_soft, ...
    CI_dpp_ens_soft(1), CI_dpp_ens_soft(2), ...
    GM_dbarpp, GSD_dbarpp, CI_dbarpp(1), CI_dbarpp(2), ...
    GM_sigmapp, GSD_sigmapp, CI_sigmapp(1), CI_sigmapp(2), ...
    GM_da, GSD_da, CI_da(1), CI_da(2), ...
    'VariableNames', {'entry_id', 'entry_label', ...
    'n_primary_particle_aggregates', 'n_morphology_aggregates', ...
    'n_primaries', 'coverage_threshold', 'logistic_bandwidth', ...
    'GM_dpp_ens_unweighted', 'GSD_dpp_ens_unweighted', ...
    'CI95_low_dpp_ens_unweighted', 'CI95_high_dpp_ens_unweighted', ...
    'GM_dpp_ens_hard', 'GSD_dpp_ens_hard', ...
    'CI95_low_dpp_ens_hard', 'CI95_high_dpp_ens_hard', ...
    'GM_dpp_ens_soft', 'GSD_dpp_ens_soft', ...
    'CI95_low_dpp_ens_soft', 'CI95_high_dpp_ens_soft', ...
    'GM_dbarpp', 'GSD_dbarpp', ...
    'CI95_low_dbarpp', 'CI95_high_dbarpp', ...
    'GM_sigmapp', 'GSD_sigmapp', ...
    'CI95_low_sigmapp', 'CI95_high_sigmapp', ...
    'GM_da', 'GSD_da', 'CI95_low_da', 'CI95_high_da'});

end

function model_inputs = build_model_inputs(entry, aggregate_table)

if ~entry.include_in_model_inputs
    model_inputs = table();
    return
end

model_table = apply_model_filter(entry, aggregate_table);
if isempty(model_table)
    warning('PFAL:main_tem_analysis_v1:EmptyModelInputs', ...
        'Entry "%s" produced no model-input rows after filtering.', entry.id);
    model_inputs = table();
    return
end

[GM_dpp, GSD_dpp, CI_dpp] = UTILS.GEOMSTATS(model_table.dbarpp_nm);
[GM_da, GSD_da, CI_da] = UTILS.GEOMSTATS(model_table.da_nm);

model_inputs = table( ...
    {entry.id}, ...
    {entry.label}, ...
    height(model_table), ...
    GM_dpp, ...
    GSD_dpp, ...
    CI_dpp(1), ...
    CI_dpp(2), ...
    GM_da, ...
    GSD_da, ...
    CI_da(1), ...
    CI_da(2), ...
    'VariableNames', {'entry_id', 'entry_label', 'n_aggregates', ...
    'GM_dpp', 'GSD_dpp', 'CI95_low_dpp', 'CI95_high_dpp', ...
    'GM_da', 'GSD_da', 'CI95_low_da', 'CI95_high_da'});

end

function model_table = apply_model_filter(entry, aggregate_table)

model_table = aggregate_table;
if isempty(fieldnames(entry.model_filter))
    return
end

field_name = entry.model_filter.field;
if ~ismember(field_name, aggregate_table.Properties.VariableNames)
    error('PFAL:main_tem_analysis_v1:InvalidModelFilterField', ...
        'Entry "%s" model_filter field "%s" is not available.', ...
        entry.id, field_name);
end

x = aggregate_table.(field_name);
v = entry.model_filter.value;
switch entry.model_filter.operator
    case '=='
        mask = x == v;
    case '~='
        mask = x ~= v;
    case '>'
        mask = x > v;
    case '>='
        mask = x >= v;
    case '<'
        mask = x < v;
    case '<='
        mask = x <= v;
    otherwise
        error('PFAL:main_tem_analysis_v1:InvalidModelFilterOperator', ...
            'Unsupported model_filter operator "%s" for entry "%s".', ...
            entry.model_filter.operator, entry.id);
end
model_table = aggregate_table(mask, :);

end

function plot_artifacts = build_plot_artifacts(entry_results, cfg_tem)

entry_styles = resolve_entry_styles(entry_results);

plot_artifacts = struct();
plot_artifacts.entry_styles = entry_styles;
plot_artifacts.ensemble_primary_summary = ...
    build_ensemble_primary_summary(entry_results);
plot_artifacts.per_aggregate_dpp_summary = build_aggregate_metric_summary( ...
    entry_results, 'dbarpp_nm', 'per_aggregate_mean_dpp');
plot_artifacts.per_aggregate_sigmapp_summary = ...
    build_aggregate_metric_summary(entry_results, 'sigmapp', ...
    'per_aggregate_sigmapp');
plot_artifacts.hybridity_frequency_summary = ...
    build_binned_frequency_summary(entry_results, 'n_subagg', ...
    cfg_tem.plots.frequency_bins.hybridity, 'subaggregate_count');
plot_artifacts.collapse_frequency_summary = ...
    build_binned_frequency_summary(entry_results, 'collapse_fraction', ...
    cfg_tem.plots.frequency_bins.collapse, 'collapse_fraction');
plot_artifacts.subaggregate_count_distribution_summary = ...
    build_exact_subaggregate_count_summary(entry_results);
plot_artifacts.plot_diagnostics = build_plot_diagnostics( ...
    entry_results, entry_styles, cfg_tem);

end

function plot_tem_analysis(entry_results, cfg_tem, plot_artifacts)

% Each figure owns an independent condition list. This design allows a
% condition to be added to one manuscript panel without changing the data
% selection or layout of any other panel.
if figure_enabled(cfg_tem, 'dpp_vs_da_manual_tem')
    plot_dpp_vs_da_publication(entry_results, cfg_tem, plot_artifacts)
end
if figure_enabled(cfg_tem, 'primary_particle_distributions')
    plot_primary_particle_distributions(entry_results, cfg_tem, ...
        plot_artifacts, 'primary_particle_distributions')
end
if figure_enabled(cfg_tem, 'appendix_a_primary_particle_distributions')
    plot_primary_particle_distributions(entry_results, cfg_tem, ...
        plot_artifacts, 'appendix_a_primary_particle_distributions')
end
if figure_enabled(cfg_tem, 'aggregate_metric_distributions')
    plot_aggregate_metric_distributions(entry_results, cfg_tem, plot_artifacts)
end
if figure_enabled(cfg_tem, 'subaggregate_count_distribution')
    plot_subaggregate_count_distribution(entry_results, cfg_tem, ...
        plot_artifacts)
end
if figure_enabled(cfg_tem, 'subaggregate_count_frequencies')
    plot_subaggregate_count_frequencies(entry_results, cfg_tem, plot_artifacts)
end
if figure_enabled(cfg_tem, 'collapsed_subaggregate_frequencies')
    plot_collapsed_subaggregate_frequencies(entry_results, cfg_tem, ...
        plot_artifacts)
end
if figure_enabled(cfg_tem, 'dpp_vs_da_by_hybridity')
    plot_dpp_vs_da_by_hybridity(entry_results, cfg_tem, plot_artifacts)
end

end

function plot_dpp_vs_da_publication(entry_results, cfg_tem, plot_artifacts)

[entry_results, styles] = select_figure_entries(entry_results, ...
    plot_artifacts.entry_styles, cfg_tem, 'dpp_vs_da_manual_tem');
assert_required_plot_metric(entry_results, 'da_nm', 'dpp_vs_da_manual_tem')
assert_required_plot_metric(entry_results, 'dbarpp_nm', 'dpp_vs_da_manual_tem')

figure_cfg = cfg_tem.plots.figures.dpp_vs_da_manual_tem;
[fig, ax, font_name] = create_publication_figure(cfg_tem, figure_cfg);

% The reference curve provides the same universal correlation used in the
% validation figures and in the PFAL scaling workflow.
[da0, dpp0] = reference_dpp_da_curve();
reference_cfg = cfg_tem.plots.reference;
plt_ref = plot(ax, da0, dpp0, 'Color', hex_to_rgb(reference_cfg.color), ...
    'LineStyle', reference_cfg.line_style, ...
    'LineWidth', reference_cfg.line_width);
hold(ax, 'on')

entry_handles = gobjects(numel(entry_results), 1);
legend_entries = cell(numel(entry_results), 1);
all_da = [];
all_dpp = [];
for i = 1:numel(entry_results)
    agg_table = entry_results(i).aggregate_table;
    style = styles(i);
    marker_size = style.scatter_size;
    entry_handles(i) = scatter(ax, agg_table.da_nm, agg_table.dbarpp_nm, ...
        marker_size, hex_to_rgb(style.scatter_color), style.marker, ...
        'MarkerFaceColor', 'none', ...
        'LineWidth', style.line_width);
    legend_entries{i} = tex_plain_text(sprintf('%s (n = %d)', ...
        entry_results(i).entry.label, height(agg_table)));
    all_da = [all_da; agg_table.da_nm(:)]; %#ok<AGROW>
    all_dpp = [all_dpp; agg_table.dbarpp_nm(:)]; %#ok<AGROW>
end

configure_axes(ax, cfg_tem, figure_cfg)
xlim(configured_axis_limit(cfg_tem, 'dpp_vs_da_manual_tem', 'xlim', ...
    nice_log_limits(all_da, [0.85 1.2]), all_da, 'log'))
ylim(configured_axis_limit(cfg_tem, 'dpp_vs_da_manual_tem', 'ylim', ...
    nice_log_limits(all_dpp, [0.95 1.05]), all_dpp, 'log'))
xticks(configured_axis_ticks(cfg_tem, 'dpp_vs_da_manual_tem', 'xticks', ...
    atems_da_ticks(), xlim(ax), 'log'))
xlabel(ax, scientific_axis_label('da', cfg_tem.plots.font.interpreter), ...
    'Interpreter', cfg_tem.plots.font.interpreter, ...
    'FontName', font_name, 'FontSize', cfg_tem.plots.font.label_size)
ylabel(ax, scientific_axis_label('dpp', cfg_tem.plots.font.interpreter), ...
    'Interpreter', cfg_tem.plots.font.interpreter, ...
    'FontName', font_name, 'FontSize', cfg_tem.plots.font.label_size)
lgd = legend(ax, [entry_handles; plt_ref], ...
    [legend_entries; {tex_plain_text('Olfert and Rogak (2019)')}], ...
    'Interpreter', figure_cfg.legend.interpreter, ...
    'FontName', resolve_font(cfg_tem.plots.font), ...
    'FontSize', cfg_tem.plots.font.legend_size, ...
    'Location', figure_cfg.legend.location, ...
    'Orientation', figure_cfg.legend.orientation, ...
    'NumColumns', figure_cfg.legend.columns, ...
    'Box', logical_to_on_off(figure_cfg.legend.box));
set(lgd, 'FontWeight', cfg_tem.plots.font.weight)

export_figure(fig, cfg_tem, figure_cfg.file_name)

end

function plot_primary_particle_distributions(entry_results, cfg_tem, ...
    plot_artifacts, figure_id)

[entry_results, styles] = select_figure_entries(entry_results, ...
    plot_artifacts.entry_styles, cfg_tem, figure_id);
assert_required_plot_metric(entry_results, 'dpp_nm', ...
    figure_id)
assert_required_plot_metric(entry_results, 'dbarpp_nm', ...
    figure_id)
assert_required_plot_metric(entry_results, 'sigmapp', ...
    figure_id)

figure_cfg = cfg_tem.plots.figures.(figure_id);
[fig, ~, ~] = create_publication_figure(cfg_tem, figure_cfg, false);
layout = tiledlayout(fig, 3, 1, 'Padding', figure_cfg.layout_padding, ...
    'TileSpacing', figure_cfg.tile_spacing);

% Panel 1 represents the condition-level primary-particle reservoir after
% the configured coverage correction. The record figure includes the KDE;
% the paper-faithful Appendix A version retains only the notched boxplots.
ax = nexttile(layout);
plot_weighted_primary_box_and_kde(entry_results, styles, cfg_tem, ...
    figure_cfg)
all_dpp = collect_primary_metric(entry_results, 'dpp_nm');
ylim(ax, configured_axis_limit(cfg_tem, figure_id, ...
    'ensemble_dpp_ylim', [5 55], all_dpp, 'log'))
yticks(ax, configured_axis_ticks(cfg_tem, figure_id, ...
    'ensemble_dpp_yticks', [5 10 20 40 80], ylim(ax), 'log'))

% Panel 2 shows the geometric mean calculated independently within each
% aggregate and therefore uses aggregate counts in the condition labels.
nexttile(layout)
    plot_publication_box_metric(entry_results, styles, 'dbarpp_nm', ...
    scientific_axis_label('dpp', cfg_tem.plots.font.interpreter), ...
    figure_cfg.aggregate_box_width, cfg_tem, 'aggregate_table', false)
all_dbarpp = collect_aggregate_metric(entry_results, 'dbarpp_nm');
ylim(configured_axis_limit(cfg_tem, figure_id, ...
    'aggregate_mean_dpp_ylim', [8 24], all_dbarpp, 'linear'))

% Panel 3 reports the within-aggregate geometric standard deviation, which
% is the direct measure of primary-particle uniformity used in the paper.
nexttile(layout)
    plot_publication_box_metric(entry_results, styles, 'sigmapp', ...
    scientific_axis_label('sigmapp', cfg_tem.plots.font.interpreter), ...
    figure_cfg.aggregate_box_width, cfg_tem, 'aggregate_table', false)
all_sigmapp = collect_aggregate_metric(entry_results, 'sigmapp');
ylim(configured_axis_limit(cfg_tem, figure_id, ...
    'sigmapp_ylim', [1.11 1.62], all_sigmapp, 'linear'))
yticks(configured_axis_ticks(cfg_tem, figure_id, ...
    'sigmapp_yticks', [1.2 1.3 1.4 1.5 1.6], ylim, 'linear'))

export_figure(fig, cfg_tem, figure_cfg.file_name)

end

function plot_aggregate_metric_distributions(entry_results, cfg_tem, ...
    plot_artifacts)

[entry_results, styles] = select_figure_entries(entry_results, ...
    plot_artifacts.entry_styles, cfg_tem, 'aggregate_metric_distributions');
assert_required_morphology_metric(entry_results, 'da_nm', ...
    'aggregate_metric_distributions')

figure_cfg = cfg_tem.plots.figures.aggregate_metric_distributions;
[fig, ~, ~] = create_publication_figure(cfg_tem, figure_cfg, false);
layout = tiledlayout(fig, 2, 2, 'Padding', figure_cfg.layout_padding, ...
    'TileSpacing', figure_cfg.tile_spacing);

% The 2x2 layout separates aggregate size from morphology. Explicitly
% nonpositive subaggregate counts are excluded; missing subaggregate counts
% do not exclude otherwise valid aggregate-scale morphology measurements.
nexttile(layout)
plot_publication_box_metric(entry_results, styles, 'da_nm', ...
    scientific_axis_label('da', cfg_tem.plots.font.interpreter), ...
    figure_cfg.box_width, cfg_tem, 'morphology_table', true)
set(gca, 'YScale', 'log')
all_da = collect_morphology_metric(entry_results, 'da_nm', true);
ylim(configured_axis_limit(cfg_tem, 'aggregate_metric_distributions', ...
    'da_ylim', nice_log_limits(all_da, [0.85 1.15]), all_da, 'log'))

nexttile(layout)
plot_publication_box_metric(entry_results, styles, 'ca', ...
    scientific_axis_label('ca', cfg_tem.plots.font.interpreter), ...
    figure_cfg.box_width, cfg_tem, 'morphology_table', true)
all_ca = collect_morphology_metric(entry_results, 'ca', true);
ylim(configured_axis_limit(cfg_tem, 'aggregate_metric_distributions', ...
    'ca_ylim', nice_linear_limits(all_ca, 0.08, [0 1]), all_ca, 'linear'))

nexttile(layout)
plot_publication_box_metric(entry_results, styles, 'zbar_opt', ...
    scientific_axis_label('zbar', cfg_tem.plots.font.interpreter), ...
    figure_cfg.box_width, cfg_tem, 'morphology_table', true)
all_zbar = collect_morphology_metric(entry_results, 'zbar_opt', true);
ylim(configured_axis_limit(cfg_tem, 'aggregate_metric_distributions', ...
    'zbar_opt_ylim', nice_linear_limits(all_zbar, 0.08, []), all_zbar, ...
    'linear'))

nexttile(layout)
plot_publication_box_metric(entry_results, styles, 'sbar_opt', ...
    scientific_axis_label('sbar', cfg_tem.plots.font.interpreter), ...
    figure_cfg.box_width, cfg_tem, 'morphology_table', true)
all_sbar = collect_morphology_metric(entry_results, 'sbar_opt', true);
ylim(configured_axis_limit(cfg_tem, 'aggregate_metric_distributions', ...
    'sbar_opt_ylim', nice_linear_limits(all_sbar, 0.08, []), all_sbar, ...
    'linear'))

export_figure(fig, cfg_tem, figure_cfg.file_name)

end

function plot_subaggregate_count_distribution(entry_results, cfg_tem, ...
    plot_artifacts)

figure_id = 'subaggregate_count_distribution';
[entry_results, styles] = select_figure_entries(entry_results, ...
    plot_artifacts.entry_styles, cfg_tem, figure_id);
assert_required_morphology_metric(entry_results, 'n_subagg', figure_id)

figure_cfg = cfg_tem.plots.figures.(figure_id);
[support, frequency_values] = exact_count_frequency_matrix( ...
    plot_artifacts.subaggregate_count_distribution_summary, entry_results);
[fig, ax, font_name] = create_publication_figure(cfg_tem, figure_cfg);
hold(ax, 'on')

legend_handles = gobjects(2, 1);
signed_frequency = [-frequency_values(:, 1), frequency_values(:, 2)];
for i = 1:2
    condition_color = hex_to_rgb(styles(i).color);
    present = frequency_values(:, i) > 0;
    present_support = support(present);
    present_frequency = signed_frequency(present, i);
    for j = 1:numel(present_support)
        line(ax, [0 present_frequency(j)], ...
            [present_support(j) present_support(j)], ...
            'Color', condition_color, ...
            'LineWidth', figure_cfg.stem_line_width, ...
            'HandleVisibility', 'off')
    end
    plot(ax, present_frequency, present_support, 'o', ...
        'LineStyle', 'none', ...
        'Color', condition_color, ...
        'MarkerSize', figure_cfg.marker_size, ...
        'MarkerFaceColor', configured_color(figure_cfg.marker_face_color), ...
        'LineWidth', figure_cfg.marker_edge_width, ...
        'HandleVisibility', 'off')
    legend_handles(i) = plot(ax, NaN, NaN, '-o', ...
        'Color', condition_color, ...
        'LineWidth', figure_cfg.stem_line_width, ...
        'MarkerSize', figure_cfg.marker_size, ...
        'MarkerFaceColor', configured_color(figure_cfg.marker_face_color));
end

configure_axes(ax, cfg_tem, figure_cfg)
if isempty(figure_cfg.x_limits)
    frequency_limit = max(10, 10 * ceil(max(frequency_values, [], 'all') / 10));
    xlim(ax, [-frequency_limit frequency_limit])
end
if isempty(figure_cfg.y_limits)
    ylim(ax, [min(support) - 0.5, max(support) + 0.5])
end
apply_configured_axis(ax, figure_cfg)
if ~isempty(figure_cfg.x_minor_ticks)
    ax.XAxis.MinorTickValues = figure_cfg.x_minor_ticks;
end
apply_mirrored_frequency_grid(ax, figure_cfg.grid)
xline(ax, 0, figure_cfg.zero_line.line_style, ...
    'Color', hex_to_rgb(figure_cfg.zero_line.color), ...
    'LineWidth', figure_cfg.zero_line.line_width, ...
    'HandleVisibility', 'off');

tick_values = ax.XTick;
tick_labels = arrayfun(@(value) sprintf('%g', abs(value)), ...
    tick_values, 'UniformOutput', false);
set(ax, 'XTickLabel', tick_labels)
xlabel(ax, 'Frequency [%]', ...
    'Interpreter', cfg_tem.plots.font.interpreter, ...
    'FontName', font_name, 'FontSize', cfg_tem.plots.font.label_size)
ylabel(ax, scientific_axis_label('n_subagg', ...
    cfg_tem.plots.font.interpreter), ...
    'Interpreter', cfg_tem.plots.font.interpreter, ...
    'FontName', font_name, 'FontSize', cfg_tem.plots.font.label_size)

legend_labels = arrayfun(@(entry) tex_plain_text(entry.entry.label), ...
    entry_results, 'UniformOutput', false);
lgd = legend(ax, legend_handles, legend_labels, ...
    'Interpreter', figure_cfg.legend.interpreter, ...
    'FontName', font_name, ...
    'FontSize', cfg_tem.plots.font.legend_size, ...
    'FontWeight', cfg_tem.plots.font.weight, ...
    'Location', figure_cfg.legend.location, ...
    'Orientation', figure_cfg.legend.orientation, ...
    'NumColumns', figure_cfg.legend.columns, ...
    'Box', logical_to_on_off(figure_cfg.legend.box));
lgd.ItemTokenSize = figure_cfg.legend_item_token_size;

export_figure(fig, cfg_tem, figure_cfg.file_name)

end

function apply_mirrored_frequency_grid(ax, grid_cfg)

set(ax, ...
    'XGrid', logical_to_on_off(grid_cfg.x_major), ...
    'XMinorGrid', logical_to_on_off(grid_cfg.x_minor), ...
    'XMinorTick', logical_to_on_off(grid_cfg.x_minor), ...
    'YGrid', logical_to_on_off(grid_cfg.y_major), ...
    'GridColor', hex_to_rgb(grid_cfg.major_color), ...
    'MinorGridColor', hex_to_rgb(grid_cfg.minor_color), ...
    'GridAlpha', grid_cfg.major_alpha, ...
    'MinorGridAlpha', grid_cfg.minor_alpha, ...
    'GridLineStyle', grid_cfg.major_line_style, ...
    'MinorGridLineStyle', grid_cfg.minor_line_style)

end

function plot_subaggregate_count_frequencies(entry_results, cfg_tem, ...
    plot_artifacts)

figure_id = 'subaggregate_count_frequencies';
[entry_results, ~] = select_figure_entries(entry_results, ...
    plot_artifacts.entry_styles, cfg_tem, figure_id);
assert_required_morphology_metric(entry_results, 'n_subagg', figure_id)

bins = cfg_tem.plots.frequency_bins.hybridity;
[frequency_values, valid_n] = frequency_matrix( ...
    plot_artifacts.hybridity_frequency_summary, entry_results, bins);
figure_cfg = cfg_tem.plots.figures.(figure_id);
[fig, ax, font_name] = create_publication_figure(cfg_tem, figure_cfg);
x_positions = frequency_x_positions(entry_results, figure_cfg);

bar_handles = bar(ax, x_positions, frequency_values, 'stacked');
style_stacked_bars(bar_handles, condition_stack_colors( ...
    entry_results, 'hybridity'), figure_cfg.bar_width, ...
    figure_cfg.bar_edge_color, figure_cfg.bar_line_width)
format_frequency_axis(ax, {entry_results.entry}, valid_n, ...
    x_positions, cfg_tem, figure_cfg)
ylabel(ax, 'Frequency [%]', 'Interpreter', ...
    cfg_tem.plots.font.interpreter, 'FontName', font_name, ...
    'FontSize', cfg_tem.plots.font.label_size)
add_frequency_legend(ax, bins, cfg_tem, figure_cfg)
add_subaggregate_arrows(ax, bar_handles, entry_results, ...
    x_positions, figure_cfg.arrows)

export_figure(fig, cfg_tem, figure_cfg.file_name)

end

function plot_collapsed_subaggregate_frequencies(entry_results, cfg_tem, ...
    plot_artifacts)

figure_id = 'collapsed_subaggregate_frequencies';
[entry_results, ~] = select_figure_entries(entry_results, ...
    plot_artifacts.entry_styles, cfg_tem, figure_id);
assert_required_morphology_metric(entry_results, 'n_subagg', figure_id)
assert_required_morphology_metric(entry_results, 'n_colaps', figure_id)

bins = cfg_tem.plots.frequency_bins.collapse;
[frequency_values, valid_n] = frequency_matrix( ...
    plot_artifacts.collapse_frequency_summary, entry_results, bins);
figure_cfg = cfg_tem.plots.figures.(figure_id);
[fig, ax, font_name] = create_publication_figure(cfg_tem, figure_cfg);
x_positions = frequency_x_positions(entry_results, figure_cfg);

bar_handles = bar(ax, x_positions, frequency_values, 'stacked');
style_stacked_bars(bar_handles, condition_stack_colors( ...
    entry_results, 'collapse'), figure_cfg.bar_width, ...
    figure_cfg.bar_edge_color, figure_cfg.bar_line_width)
format_frequency_axis(ax, {entry_results.entry}, valid_n, ...
    x_positions, cfg_tem, figure_cfg)
ylabel(ax, 'Frequency [%]', 'Interpreter', ...
    cfg_tem.plots.font.interpreter, 'FontName', font_name, ...
    'FontSize', cfg_tem.plots.font.label_size)
add_frequency_legend(ax, bins, cfg_tem, figure_cfg)

export_figure(fig, cfg_tem, figure_cfg.file_name)

end

function x_positions = frequency_x_positions(entry_results, figure_cfg)

x_positions = 1 + (0:numel(entry_results) - 1) * ...
    figure_cfg.condition_spacing;

end

function add_frequency_legend(ax, bins, cfg_tem, figure_cfg)

lgd = legend(ax, {bins.label}, ...
    'Interpreter', figure_cfg.legend.interpreter, ...
    'FontName', resolve_font(cfg_tem.plots.font), ...
    'FontSize', cfg_tem.plots.font.legend_size, ...
    'FontWeight', cfg_tem.plots.font.weight, ...
    'Location', figure_cfg.legend.location, ...
    'Orientation', figure_cfg.legend.orientation, ...
    'NumColumns', figure_cfg.legend.columns, ...
    'Box', logical_to_on_off(figure_cfg.legend.box));
lgd.ItemTokenSize = figure_cfg.legend_item_token_size;

end

function plot_dpp_vs_da_by_hybridity(entry_results, cfg_tem, plot_artifacts)

[entry_results, ~] = select_figure_entries(entry_results, ...
    plot_artifacts.entry_styles, cfg_tem, 'dpp_vs_da_by_hybridity');
assert_required_plot_metric(entry_results, 'n_subagg', ...
    'dpp_vs_da_by_hybridity')
assert_required_plot_metric(entry_results, 'da_nm', ...
    'dpp_vs_da_by_hybridity')
assert_required_plot_metric(entry_results, 'dbarpp_nm', ...
    'dpp_vs_da_by_hybridity')

figure_cfg = cfg_tem.plots.figures.dpp_vs_da_by_hybridity;
[fig, ax, font_name] = create_publication_figure(cfg_tem, figure_cfg);

[da0, dpp0] = reference_dpp_da_curve();
reference_cfg = cfg_tem.plots.reference;
plt_ref = plot(ax, da0, dpp0, 'Color', hex_to_rgb(reference_cfg.color), ...
    'LineStyle', reference_cfg.line_style, ...
    'LineWidth', reference_cfg.line_width);
hold(ax, 'on')

all_da = collect_aggregate_metric(entry_results, 'da_nm');
all_dpp = collect_aggregate_metric(entry_results, 'dbarpp_nm');
all_n_hyb = collect_aggregate_metric(entry_results, 'n_subagg');

categories = figure_cfg.categories;
category_handles = gobjects(numel(categories), 1);
for i = 1:numel(categories)
    mask = interval_mask(all_n_hyb, categories(i));
    category_handles(i) = scatter(ax, all_da(mask), all_dpp(mask), ...
        categories(i).marker_size, hex_to_rgb(categories(i).color), ...
        categories(i).marker, 'MarkerFaceColor', 'none', ...
        'LineWidth', categories(i).line_width);
end

configure_axes(ax, cfg_tem, figure_cfg)
xlim(configured_axis_limit(cfg_tem, 'dpp_vs_da_by_hybridity', 'xlim', ...
    nice_log_limits(all_da, [0.8 1.2]), all_da, 'log'))
ylim(configured_axis_limit(cfg_tem, 'dpp_vs_da_by_hybridity', 'ylim', ...
    nice_log_limits(all_dpp, [0.95 1.05]), all_dpp, 'log'))
xlabel(ax, scientific_axis_label('da', cfg_tem.plots.font.interpreter), ...
    'Interpreter', cfg_tem.plots.font.interpreter, 'FontName', font_name, ...
    'FontSize', cfg_tem.plots.font.label_size)
ylabel(ax, scientific_axis_label('dpp', cfg_tem.plots.font.interpreter), ...
    'Interpreter', cfg_tem.plots.font.interpreter, 'FontName', font_name, ...
    'FontSize', cfg_tem.plots.font.label_size)
lgd = legend(ax, [plt_ref; category_handles], ...
    [{tex_plain_text('Olfert and Rogak (2019)')}, {categories.label}], ...
    'Interpreter', figure_cfg.legend.interpreter, ...
    'FontName', resolve_font(cfg_tem.plots.font), ...
    'FontSize', cfg_tem.plots.font.legend_size, ...
    'Location', figure_cfg.legend.location, ...
    'Orientation', figure_cfg.legend.orientation, ...
    'NumColumns', figure_cfg.legend.columns, ...
    'Box', logical_to_on_off(figure_cfg.legend.box));
set(lgd, 'FontWeight', cfg_tem.plots.font.weight)

export_figure(fig, cfg_tem, figure_cfg.file_name)

end

function plot_weighted_primary_box_and_kde(entry_results, styles, cfg_tem, ...
    figure_cfg)

values = [];
groups = [];
labels = condition_labels(entry_results, 'primary');
n_resample = cfg_tem.analysis.random_resample_count;

for i = 1:numel(entry_results)
    primary_table = entry_results(i).primary_table;
    x = primary_table.dpp_nm;
    w = primary_table.weight_soft;
    mask = isfinite(x) & x > 0 & isfinite(w) & w > 0;
    x = x(mask);
    w = w(mask);
    if isempty(x)
        continue
    end
    idx = weighted_sample_indices(numel(x), n_resample, w);
    values = [values; x(idx)]; %#ok<AGROW>
    groups = [groups; repmat(i, n_resample, 1)]; %#ok<AGROW>
end

if isempty(values)
    error('PFAL:main_tem_analysis_v1:MissingPrimaryDistribution', ...
        'No weighted primary-particle values are available for plotting.')
end

boxplot(values, groups, 'Labels', labels, ...
    'Notch', logical_to_on_off(cfg_tem.plots.boxplot.notch), ...
    'Symbol', cfg_tem.plots.boxplot.symbol, ...
    'Widths', figure_cfg.ensemble_box_width, ...
    'Colors', box_line_colors(styles))
style_boxplot(gca, styles, cfg_tem.plots.boxplot)
hold on

if figure_cfg.kde_enabled
    for i = 1:numel(entry_results)
        primary_table = entry_results(i).primary_table;
        x = primary_table.dpp_nm;
        w = primary_table.weight_soft;
        mask = isfinite(x) & x > 0 & isfinite(w) & w > 0;
        x = x(mask);
        w = w(mask);
        if numel(x) < 2
            continue
        end
        [f, xi] = ksdensity(log10(x), 'Weights', w / sum(w));
        finite_mask = isfinite(f) & isfinite(xi);
        f = f(finite_mask);
        xi = xi(finite_mask);
        if isempty(f) || max(f) <= 0
            continue
        end
        y = 10 .^ xi;
        x_anchor = i - figure_cfg.kde_offset;
        density_x = x_anchor - figure_cfg.kde_width * f / max(f);
        plot(density_x, y, ...
            'Color', hex_to_rgb(styles(i).box_edge_color), ...
            'LineWidth', styles(i).line_width)
        fill([density_x, x_anchor * ones(size(density_x))], ...
            [y, fliplr(y)], hex_to_rgb(styles(i).box_face_color), ...
            'FaceAlpha', styles(i).box_face_alpha, ...
            'EdgeColor', 'none');
    end
end

ax = gca;
configure_axes(ax, cfg_tem, figure_cfg)
set(ax, 'YScale', 'log')
ylabel(ax, scientific_axis_label('dpp_individual', ...
    cfg_tem.plots.font.interpreter), ...
    'Interpreter', cfg_tem.plots.font.interpreter, ...
    'FontSize', cfg_tem.plots.font.label_size)
if figure_cfg.kde_enabled
    xlim([figure_cfg.group_padding, ...
        numel(entry_results) + figure_cfg.group_padding])
else
    xlim([1 - figure_cfg.group_padding, ...
        numel(entry_results) + figure_cfg.group_padding])
end

end

function plot_publication_box_metric(entry_results, styles, metric_name, ...
    y_label, box_width, cfg_tem, table_field, filter_hybrid_rows)

values = [];
groups = [];
labels = condition_labels(entry_results, table_field);
for i = 1:numel(entry_results)
    source_table = entry_results(i).(table_field);
    metric_values = source_table.(metric_name);
    if filter_hybrid_rows && ismember('n_subagg', ...
            source_table.Properties.VariableNames)
        hybrid_mask = source_table.n_subagg > 0 | ...
            ~isfinite(source_table.n_subagg);
        metric_values = metric_values(hybrid_mask);
    end
    metric_values = metric_values(isfinite(metric_values));
    values = [values; metric_values(:)]; %#ok<AGROW>
    groups = [groups; repmat(i, numel(metric_values), 1)]; %#ok<AGROW>
end

if isempty(values)
    text(0.5, 0.5, sprintf('No finite %s values', metric_name), ...
        'HorizontalAlignment', 'center', 'interpreter', 'none')
    axis off
    return
end

boxplot(values, groups, 'Labels', labels, ...
    'Notch', logical_to_on_off(cfg_tem.plots.boxplot.notch), ...
    'Symbol', cfg_tem.plots.boxplot.symbol, 'Widths', box_width, ...
    'Colors', box_line_colors(styles))
style_boxplot(gca, styles, cfg_tem.plots.boxplot)
ax = gca;
font_name = resolve_font(cfg_tem.plots.font);
set(ax, 'TickLabelInterpreter', 'tex', ...
    'FontName', font_name, 'FontSize', cfg_tem.plots.font.axis_size, ...
    'FontWeight', cfg_tem.plots.font.weight, ...
    'TickLength', cfg_tem.plots.axes.tick_length, ...
    'LineWidth', cfg_tem.plots.axis_line_width, ...
    'Box', logical_to_on_off(cfg_tem.plots.axes.box), ...
    'Layer', cfg_tem.plots.axes.layer)
ylabel(ax, y_label, 'Interpreter', cfg_tem.plots.font.interpreter, ...
    'FontName', font_name, 'FontSize', cfg_tem.plots.font.label_size)

end

function [entry_results_out, styles_out, mask] = select_figure_entries( ...
    entry_results, styles, cfg_tem, figure_id)

% Preserve the configured order because it controls left-to-right condition
% placement and legend order in the exported figure.
condition_ids = cfg_tem.plots.figures.(figure_id).condition_ids;
available_ids = arrayfun(@(x) x.entry.id, entry_results, ...
    'UniformOutput', false);
mask = false(numel(entry_results), 1);
ordered_indices = zeros(numel(condition_ids), 1);
for i = 1:numel(condition_ids)
    index = find(strcmp(available_ids, condition_ids{i}), 1);
    if isempty(index)
        error('PFAL:main_tem_analysis_v1:ConditionNotProcessed', ...
            ['Figure "%s" requests condition "%s", but that condition ', ...
            'was not processed. Confirm that the entry is enabled.'], ...
            figure_id, condition_ids{i});
    end
    ordered_indices(i) = index;
    mask(index) = true;
end
entry_results_out = entry_results(ordered_indices);
styles_out = styles(ordered_indices);

end

function assert_required_plot_metric(entry_results, metric_name, figure_name)

has_values = false;
for i = 1:numel(entry_results)
    if ismember(metric_name, entry_results(i).aggregate_table.Properties.VariableNames)
        values = entry_results(i).aggregate_table.(metric_name);
    elseif ismember(metric_name, entry_results(i).primary_table.Properties.VariableNames)
        values = entry_results(i).primary_table.(metric_name);
    else
        error('PFAL:main_tem_analysis_v1:MissingPlotMetric', ...
            'Figure "%s" requires metric "%s", but entry "%s" does not contain it.', ...
            figure_name, metric_name, entry_results(i).entry.id);
    end
    has_values = has_values || any(isfinite(values));
end

if ~has_values
    error('PFAL:main_tem_analysis_v1:EmptyPlotMetric', ...
        'Figure "%s" requires metric "%s", but no finite values are available.', ...
        figure_name, metric_name);
end

end

function assert_required_morphology_metric(entry_results, metric_name, ...
    figure_name)

has_values = false;
for i = 1:numel(entry_results)
    morphology_table = entry_results(i).morphology_table;
    if ~ismember(metric_name, morphology_table.Properties.VariableNames)
        error('PFAL:main_tem_analysis_v1:MissingMorphologyMetric', ...
            ['Figure "%s" requires morphology metric "%s", but entry ', ...
            '"%s" does not contain it.'], ...
            figure_name, metric_name, entry_results(i).entry.id);
    end
    has_values = has_values || any(isfinite(morphology_table.(metric_name)));
end
if ~has_values
    error('PFAL:main_tem_analysis_v1:EmptyMorphologyMetric', ...
        ['Figure "%s" requires morphology metric "%s", but no finite ', ...
        'values are available.'], figure_name, metric_name);
end

end

function styles = resolve_entry_styles(entry_results)

styles = repmat(default_style(), numel(entry_results), 1);
for i = 1:numel(entry_results)
    entry = entry_results(i).entry;
    source = entry.style;
    style = default_style();
    style.entry_id = entry.id;
    style.entry_label = entry.label;
    style.entry_type = entry.entry_type;
    style.color = source.color;
    style.secondary_color = source.secondary_color;
    style.scatter_color = source.color;
    style.box_face_color = entry.style.box_face_color;
    style.box_edge_color = source.secondary_color;
    style.median_color = source.median_color;
    style.marker = source.marker;
    style.scatter_size = source.marker_size;
    style.line_width = source.line_width;
    style.box_face_alpha = source.box_face_alpha;
    styles(i) = style;
end

end

function style = default_style()

style = struct( ...
    'entry_id', '', ...
    'entry_label', '', ...
    'entry_type', '', ...
    'color', '', ...
    'secondary_color', '', ...
    'scatter_color', '', ...
    'box_face_color', '', ...
    'box_edge_color', '', ...
    'median_color', '', ...
    'marker', 'o', ...
    'scatter_size', 30, ...
    'line_width', 1.5, ...
    'box_face_alpha', 0.25);

end

function style_boxplot(ax, styles, boxplot_cfg)

boxes = findobj(ax, 'Tag', 'Box');
for i = 1:numel(boxes)
    group_index = nearest_group_index(boxes(i), numel(styles));
    patch(get(boxes(i), 'XData'), get(boxes(i), 'YData'), ...
        hex_to_rgb(styles(group_index).box_face_color), ...
        'EdgeColor', hex_to_rgb(styles(group_index).box_edge_color), ...
        'FaceAlpha', styles(group_index).box_face_alpha, ...
        'LineWidth', styles(group_index).line_width);
end

medians = findobj(ax, 'Tag', 'Median');
for i = 1:numel(medians)
    group_index = nearest_group_index(medians(i), numel(styles));
    set(medians(i), 'Color', hex_to_rgb(styles(group_index).median_color), ...
        'LineWidth', max(boxplot_cfg.median_min_line_width, ...
        styles(group_index).line_width));
end

outliers = findobj(ax, 'Tag', 'Outliers');
for i = 1:numel(outliers)
    group_index = nearest_group_index(outliers(i), numel(styles));
    edge_color = hex_to_rgb(styles(group_index).box_edge_color);
    outliers(i).Color = edge_color;
    outliers(i).MarkerEdgeColor = edge_color;
    if boxplot_cfg.outlier_filled
        outliers(i).MarkerFaceColor = edge_color;
    else
        outliers(i).MarkerFaceColor = 'none';
    end
    outliers(i).MarkerSize = boxplot_cfg.outlier_marker_size;
end

set(findobj(ax, 'type', 'line', 'tag', 'Upper Whisker'), ...
    'linestyle', boxplot_cfg.whisker_line_style)
set(findobj(ax, 'type', 'line', 'tag', 'Lower Whisker'), ...
    'linestyle', boxplot_cfg.whisker_line_style)

end

function colors = box_line_colors(styles)

colors = zeros(numel(styles), 3);
for i = 1:numel(styles)
    colors(i, :) = hex_to_rgb(styles(i).box_edge_color);
end

end

function group_index = nearest_group_index(h, n_groups)

x_data = get(h, 'XData');
x_data = x_data(isfinite(x_data));
if isempty(x_data)
    group_index = 1;
else
    group_index = round(median(x_data));
end
group_index = min(max(group_index, 1), n_groups);

end

function summary_table = build_ensemble_primary_summary(entry_results)

entry_id = {};
entry_label = {};
stat_mode = {};
n_values = [];
gm_value = [];
gsd_value = [];
ci_low = [];
ci_high = [];

for i = 1:numel(entry_results)
    primary_table = entry_results(i).primary_table;
    modes = {
        'unweighted', ones(height(primary_table), 1);
        'hard_coverage_weighted', primary_table.weight_hard;
        'soft_coverage_weighted', primary_table.weight_soft};
    for j = 1:size(modes, 1)
        [gm, gsd, ci] = UTILS.GEOMSTATS(primary_table.dpp_nm, modes{j, 2});
        entry_id{end + 1, 1} = entry_results(i).entry.id; %#ok<AGROW>
        entry_label{end + 1, 1} = entry_results(i).entry.label; %#ok<AGROW>
        stat_mode{end + 1, 1} = modes{j, 1}; %#ok<AGROW>
        n_values(end + 1, 1) = height(primary_table); %#ok<AGROW>
        gm_value(end + 1, 1) = gm; %#ok<AGROW>
        gsd_value(end + 1, 1) = gsd; %#ok<AGROW>
        ci_low(end + 1, 1) = ci(1); %#ok<AGROW>
        ci_high(end + 1, 1) = ci(2); %#ok<AGROW>
    end
end

summary_table = table(entry_id, entry_label, stat_mode, n_values, ...
    gm_value, gsd_value, ci_low, ci_high, 'VariableNames', ...
    {'entry_id', 'entry_label', 'stat_mode', 'n_primaries', ...
    'GM_dpp', 'GSD_dpp', 'CI95_low_dpp', 'CI95_high_dpp'});

end

function summary_table = build_aggregate_metric_summary(entry_results, ...
    metric_name, metric_label)

entry_id = {};
entry_label = {};
n_values = [];
gm_value = [];
gsd_value = [];
ci_low = [];
ci_high = [];

for i = 1:numel(entry_results)
    values = entry_results(i).aggregate_table.(metric_name);
    values = values(isfinite(values) & values > 0);
    [gm, gsd, ci] = UTILS.GEOMSTATS(values);
    entry_id{end + 1, 1} = entry_results(i).entry.id; %#ok<AGROW>
    entry_label{end + 1, 1} = entry_results(i).entry.label; %#ok<AGROW>
    n_values(end + 1, 1) = numel(values); %#ok<AGROW>
    gm_value(end + 1, 1) = gm; %#ok<AGROW>
    gsd_value(end + 1, 1) = gsd; %#ok<AGROW>
    ci_low(end + 1, 1) = ci(1); %#ok<AGROW>
    ci_high(end + 1, 1) = ci(2); %#ok<AGROW>
end

summary_table = table(entry_id, entry_label, repmat({metric_label}, ...
    numel(entry_id), 1), n_values, gm_value, gsd_value, ci_low, ...
    ci_high, 'VariableNames', {'entry_id', 'entry_label', ...
    'metric', 'n_aggregates', 'GM_value', 'GSD_value', ...
    'CI95_low_value', 'CI95_high_value'});

end

function summary_table = build_exact_subaggregate_count_summary(entry_results)

% Include zero-frequency integer counts so the exported table reproduces
% the complete support used by the mirrored distribution.
valid_values = cell(numel(entry_results), 1);
maximum_count = 0;
for i = 1:numel(entry_results)
    values = double(entry_results(i).morphology_table.n_subagg);
    values = values(isfinite(values) & values >= 1 & values == round(values));
    valid_values{i} = values(:);
    if ~isempty(values)
        maximum_count = max(maximum_count, max(values));
    end
end
if maximum_count < 1
    summary_table = table(cell(0, 1), cell(0, 1), zeros(0, 1), ...
        zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
        'VariableNames', {'entry_id','entry_label','n_subagg', ...
        'valid_n','count','frequency_percent'});
    return
end

support = (1:maximum_count).';
rows = cell(numel(entry_results), 1);
for i = 1:numel(entry_results)
    values = valid_values{i};
    counts = arrayfun(@(value) nnz(values == value), support);
    valid_n = numel(values);
    rows{i} = table( ...
        repmat({entry_results(i).entry.id}, maximum_count, 1), ...
        repmat({entry_results(i).entry.label}, maximum_count, 1), ...
        support, repmat(valid_n, maximum_count, 1), counts, ...
        100 * counts / max(valid_n, 1), ...
        'VariableNames', {'entry_id','entry_label','n_subagg', ...
        'valid_n','count','frequency_percent'});
end
summary_table = vertcat(rows{:});

end

function frequency_table = build_binned_frequency_summary(entry_results, ...
    value_source, bins, metric_name)

% Store one row per condition and bin. This long-form schema remains stable
% when bins are added, removed, or renamed in a configuration file.
rows = cell(numel(entry_results) * numel(bins), 1);
row_index = 0;
for i = 1:numel(entry_results)
    values = frequency_metric_values(entry_results(i), value_source);
    assignments = zeros(size(values));
    for j = 1:numel(bins)
        mask = interval_mask(values, bins(j));
        assignments = assignments + mask;
        row_index = row_index + 1;
        rows{row_index} = table( ...
            {entry_results(i).entry.id}, {entry_results(i).entry.label}, ...
            {metric_name}, {bins(j).id}, {bins(j).label}, ...
            bins(j).lower, bins(j).upper, bins(j).include_lower, ...
            bins(j).include_upper, numel(values), nnz(mask), ...
            100 * nnz(mask) / max(numel(values), 1), ...
            'VariableNames', {'entry_id','entry_label','metric', ...
            'bin_id','bin_label','lower_bound','upper_bound', ...
            'include_lower','include_upper','valid_n','count', ...
            'frequency_percent'});
    end
    if any(assignments ~= 1)
        error('PFAL:main_tem_analysis_v1:IncompleteFrequencyBins', ...
            ['Configured bins for metric "%s" must assign every valid ', ...
            'value exactly once.'], metric_name);
    end
end
frequency_table = vertcat(rows{1:row_index});

end

function values = frequency_metric_values(entry_result, value_source)

switch value_source
    case 'n_subagg'
        values = entry_result.morphology_table.n_subagg;
        values = values(isfinite(values) & values >= 1);
    case 'collapse_fraction'
        values = collapse_fraction_values(entry_result.morphology_table, ...
            entry_result.entry.id);
    otherwise
        error('PFAL:main_tem_analysis_v1:UnknownFrequencyMetric', ...
            'Unknown frequency metric source "%s".', value_source);
end

end

function diagnostics = build_plot_diagnostics(entry_results, styles, cfg_tem)

entry_id = {};
entry_label = {};
entry_type = {};
scatter_color = {};
box_face_color = {};
box_edge_color = {};
median_color = {};
marker = {};
missing_optional_metrics = {};
diagnostic = {};
optional_metrics = {'n_subagg', 'n_colaps', 'ca', 'zbar_opt', 'sbar_opt'};

for i = 1:numel(entry_results)
    agg_table = entry_results(i).morphology_table;
    missing = {};
    for j = 1:numel(optional_metrics)
        metric = optional_metrics{j};
        if ~ismember(metric, agg_table.Properties.VariableNames) || ...
                ~any(isfinite(agg_table.(metric)))
            missing{end + 1} = metric; %#ok<AGROW>
        end
    end

    entry_id{end + 1, 1} = entry_results(i).entry.id; %#ok<AGROW>
    entry_label{end + 1, 1} = entry_results(i).entry.label; %#ok<AGROW>
    entry_type{end + 1, 1} = entry_results(i).entry.entry_type; %#ok<AGROW>
    scatter_color{end + 1, 1} = styles(i).scatter_color; %#ok<AGROW>
    box_face_color{end + 1, 1} = styles(i).box_face_color; %#ok<AGROW>
    box_edge_color{end + 1, 1} = styles(i).box_edge_color; %#ok<AGROW>
    median_color{end + 1, 1} = styles(i).median_color; %#ok<AGROW>
    marker{end + 1, 1} = styles(i).marker; %#ok<AGROW>
    missing_optional_metrics{end + 1, 1} = strjoin(missing, ','); %#ok<AGROW>
    diagnostic{end + 1, 1} = sprintf( ...
        'Configured style; config=%s; PNG resolution=%d dpi', ...
        cfg_tem.config_file, cfg_tem.plots.png_resolution); %#ok<AGROW>
end

diagnostics = table(entry_id, entry_label, entry_type, scatter_color, ...
    box_face_color, box_edge_color, median_color, marker, ...
    missing_optional_metrics, ...
    diagnostic, 'VariableNames', {'entry_id', 'entry_label', ...
    'entry_type', 'scatter_color', 'box_face_color', 'box_edge_color', ...
    'median_color', 'marker', 'missing_optional_metrics', ...
    'diagnostic'});

end

function labels = condition_labels(entry_results, count_kind, count_override)

labels = cell(numel(entry_results), 1);
for i = 1:numel(entry_results)
    if nargin >= 3 && ~isempty(count_override)
        n_value = count_override(i);
    else
        switch count_kind
            case 'primary'
                n_value = height(entry_results(i).primary_table);
            case 'morphology_table'
                n_value = height(entry_results(i).morphology_table);
            otherwise
                n_value = height(entry_results(i).aggregate_table);
        end
    end
    labels{i} = sprintf('%s\\newline(n = %d)', ...
        latex_text(entry_results(i).entry.label), n_value);
end

end

function values = collect_aggregate_metric(entry_results, metric_name, varargin)

filter_hybrid_rows = false;
if nargin > 2
    filter_hybrid_rows = varargin{1};
end

values = [];
for i = 1:numel(entry_results)
    metric_values = entry_results(i).aggregate_table.(metric_name);
    if filter_hybrid_rows && ismember('n_subagg', ...
            entry_results(i).aggregate_table.Properties.VariableNames)
        hybrid_mask = entry_results(i).aggregate_table.n_subagg > 0 | ...
            ~isfinite(entry_results(i).aggregate_table.n_subagg);
        metric_values = metric_values(hybrid_mask);
    end
    values = [values; metric_values(:)]; %#ok<AGROW>
end
values = values(isfinite(values));

end

function values = collect_morphology_metric(entry_results, metric_name, varargin)

filter_hybrid_rows = false;
if nargin > 2
    filter_hybrid_rows = varargin{1};
end

values = [];
for i = 1:numel(entry_results)
    morphology_table = entry_results(i).morphology_table;
    metric_values = morphology_table.(metric_name);
    if filter_hybrid_rows && ismember('n_subagg', ...
            morphology_table.Properties.VariableNames)
        valid_rows = morphology_table.n_subagg > 0 | ...
            ~isfinite(morphology_table.n_subagg);
        metric_values = metric_values(valid_rows);
    end
    values = [values; metric_values(:)]; %#ok<AGROW>
end
values = values(isfinite(values));

end

function values = collect_primary_metric(entry_results, metric_name)

values = [];
for i = 1:numel(entry_results)
    values = [values; entry_results(i).primary_table.(metric_name)(:)]; %#ok<AGROW>
end
values = values(isfinite(values));

end

function [da0, dpp0] = reference_dpp_da_curve()

r0 = (2e4 / 1e0) ^ (1 / (1e4 - 1));
da0 = 1e0 * ones(1e4, 1) .* r0 .^ (((1:1e4) - 1)');
dpp0 = 17.8 * (da0 / 100) .^ 0.35;

end

function ticks = atems_da_ticks()

% The selected log-scale ticks retain detail below 100 nm without causing
% label overlap when the configured condition count changes the axes width.
ticks = [30 50 70 100 200 300 500 700];

end

function style_stacked_bars(bar_handles, colors, bar_width, ...
    edge_color, line_width)

% Apply one condition-specific color sequence to each stacked bar while
% keeping the category order consistent across conditions. The edge style
% is configured separately so borderless manuscript bars remain explicit.

if size(colors, 2) ~= numel(bar_handles)
    error('PFAL:main_tem_analysis_v1:InsufficientPaletteColors', ...
        'The configured stack palette must contain one color per bin.');
end
for i = 1:numel(bar_handles)
    bar_handles(i).BarWidth = bar_width;
    bar_handles(i).FaceColor = 'flat';
    bar_handles(i).CData = squeeze(colors(:, i, :));
    bar_handles(i).EdgeColor = configured_color(edge_color);
    if line_width > 0
        bar_handles(i).LineWidth = line_width;
    end
end

end

function color = configured_color(value)

% Preserve MATLAB's named color controls, such as "none", while converting
% configured hexadecimal colors to numeric RGB triplets.
if strcmpi(value, 'none')
    color = 'none';
else
    color = hex_to_rgb(value);
end

end

function colors = condition_stack_colors(entry_results, metric_name)

% The first dimension follows condition order, the second follows category
% order, and the third contains RGB components for MATLAB bar CData.
n_conditions = numel(entry_results);
n_bins = numel(entry_results(1).entry.style.stack_colors.(metric_name));
colors = zeros(n_conditions, n_bins, 3);
for i = 1:n_conditions
    palette = entry_results(i).entry.style.stack_colors.(metric_name);
    if numel(palette) ~= n_bins
        error('PFAL:main_tem_analysis_v1:InconsistentStackPalette', ...
            'Condition stack palettes must contain the same number of colors.');
    end
    colors(i, :, :) = reshape(color_list_to_rgb(palette), 1, n_bins, 3);
end

end

function format_frequency_axis(ax, entries, valid_n, x_positions, ...
    cfg_tem, figure_cfg)

labels = cell(numel(entries), 1);
for i = 1:numel(entries)
    labels{i} = sprintf('%s\\newline(n = %d)', ...
        latex_text(entries{i}.label), valid_n(i));
end
configure_axes(ax, cfg_tem, figure_cfg)
set(ax, 'XTick', x_positions, 'XTickLabel', labels, ...
    'TickLabelInterpreter', 'tex')
xtickangle(ax, 0)
if numel(x_positions) > 1
    side_margin = figure_cfg.condition_spacing * 0.55;
else
    side_margin = 0.5;
end
xlim(ax, [x_positions(1) - side_margin, ...
    x_positions(end) + side_margin])
apply_configured_axis(ax, figure_cfg)

end

function add_subaggregate_arrows(ax, bar_handles, entry_results, ...
    x_positions, arrow_cfg)

if isempty(bar_handles) || ~arrow_cfg.enabled
    return
end
drawnow
entry_ids = arrayfun(@(r) string(r.entry.id), entry_results);
for arrow_id = string(arrow_cfg.condition_ids(:)).'
    i = find(entry_ids == arrow_id, 1);
    if isempty(i)
        continue
    end
    target_index = arrow_cfg.target_bin_index;
    segment_values = arrayfun(@(h) h.YData(i), bar_handles);
    segment_height = segment_values(target_index);
    if ~isfinite(segment_height) || segment_height <= 0
        continue
    end
    segment_bottom = sum(segment_values(1:target_index - 1));
    % Anchor the shaft to the selected segment in data coordinates so legend
    % placement and axes resizing cannot move the arrow away from the bar.
    segment_boundary = segment_bottom + ...
        arrow_cfg.tail_fraction * segment_height;
    tail_y = min(arrow_cfg.maximum_tail_y, segment_boundary);
    head_y = max(segment_bottom + arrow_cfg.minimum_head_clearance, ...
        tail_y - arrow_cfg.shaft_length);
    if head_y >= tail_y
        continue
    end
    add_down_arrow(ax, x_positions(i), tail_y, head_y, arrow_cfg)
end

end

function add_down_arrow(ax, x_value, y_start, y_end, arrow_cfg)

arrow_color = hex_to_rgb(arrow_cfg.color);
line(ax, [x_value x_value], [y_start y_end], ...
    'Color', arrow_color, 'LineWidth', arrow_cfg.line_width, ...
    'Clipping', 'on', 'HandleVisibility', 'off')
line(ax, x_value, y_end, 'LineStyle', 'none', ...
    'Marker', arrow_cfg.head_marker, ...
    'MarkerSize', arrow_cfg.head_marker_size, ...
    'MarkerFaceColor', arrow_color, 'MarkerEdgeColor', arrow_color, ...
    'Clipping', 'on', 'HandleVisibility', 'off')

end
function tf = figure_enabled(cfg_tem, figure_id)

tf = cfg_tem.plots.enabled && cfg_tem.plots.figures.(figure_id).enabled;

end

function [fig, ax, font_name] = create_publication_figure( ...
    cfg_tem, figure_cfg, create_axes)

if nargin < 3
    create_axes = true;
end
font_name = resolve_font(cfg_tem.plots.font);
figure_visibility = cfg_tem.plots.visible;
if cfg_tem.plots.export
    % Off-screen construction prevents the monitor work area from constraining
    % the configured publication aspect ratio. The requested window visibility
    % is restored after export.
    figure_visibility = 'off';
end
fig = figure('Visible', figure_visibility, 'Color', 'white', ...
    'Position', figure_cfg.position, ...
    'DefaultAxesFontName', font_name, ...
    'DefaultTextFontName', font_name);
if create_axes
    ax = axes(fig);
else
    ax = gobjects(0);
end

end

function configure_axes(ax, cfg_tem, figure_cfg)

font_name = resolve_font(cfg_tem.plots.font);
set(ax, 'Box', logical_to_on_off(cfg_tem.plots.axes.box), ...
    'Layer', cfg_tem.plots.axes.layer, ...
    'TickLength', cfg_tem.plots.axes.tick_length, ...
    'XScale', figure_cfg.x_scale, 'YScale', figure_cfg.y_scale, ...
    'FontName', font_name, 'FontSize', cfg_tem.plots.font.axis_size, ...
    'FontWeight', cfg_tem.plots.font.weight, ...
    'TickLabelInterpreter', cfg_tem.plots.font.interpreter, ...
    'LineWidth', cfg_tem.plots.axis_line_width)
if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
    ax.Toolbar.Visible = 'off';
end

end

function apply_configured_axis(ax, figure_cfg)

if ~isempty(figure_cfg.x_limits)
    xlim(ax, figure_cfg.x_limits)
end
if ~isempty(figure_cfg.y_limits)
    ylim(ax, figure_cfg.y_limits)
end
if ~isempty(figure_cfg.x_ticks)
    xticks(ax, figure_cfg.x_ticks)
end
if ~isempty(figure_cfg.y_ticks)
    yticks(ax, figure_cfg.y_ticks)
end

end

function font_name = resolve_font(font_cfg)

available_fonts = listfonts;
match = find(strcmpi(available_fonts, font_cfg.family), 1);
if ~isempty(match)
    font_name = available_fonts{match};
    return
end

match = find(strcmpi(available_fonts, font_cfg.fallback), 1);
if ~isempty(match)
    font_name = available_fonts{match};
    warning('PFAL:main_tem_analysis_v1:FontFallback', ...
        'Font "%s" is unavailable; using "%s".', ...
        font_cfg.family, font_name);
    return
end
font_name = get(groot, 'DefaultAxesFontName');
warning('PFAL:main_tem_analysis_v1:FontFallback', ...
    'Configured fonts are unavailable; using "%s".', font_name);

end

function label = scientific_axis_label(quantity, interpreter)

switch lower(interpreter)
    case 'latex'
        labels = struct('da', '$d_{\mathrm{a}}$ [nm]', ...
            'dpp', '$d_{\mathrm{pp}}$ [nm]', ...
            'dpp_individual', '$d_{\mathrm{pp}}^{(i,j)}$ [nm]', ...
            'sigmapp', '$\sigma_{\mathrm{pp}}$ [-]', ...
            'ca', '$c_{\mathrm{a}}$ [-]', ...
            'zbar', '$z_{\mathrm{a}}$ [-]', ...
            'sbar', '$s_{\mathrm{a}}$ [-]', ...
            'n_subagg', '$N_{\mathrm{subagg}}$ [-]', ...
            'f_col_subagg', '$f_{\mathrm{col,subagg}}$ [-]');
    case 'tex'
        labels = struct('da', '{\it d}_{\rm a} [nm]', ...
            'dpp', '{\it d}_{\rm pp} [nm]', ...
            'dpp_individual', '{\it d}_{\rm pp}^{(i,j)} [nm]', ...
            'sigmapp', '{\it \sigma}_{\rm pp} [-]', ...
            'ca', '{\it c}_{\rm a} [-]', ...
            'zbar', '{\it z}_{\rm a} [-]', ...
            'sbar', '{\it s}_{\rm a} [-]', ...
            'n_subagg', '{\it N}_{\rm subagg} [-]', ...
            'f_col_subagg', '{\it f}_{\rm col,subagg} [-]');
    otherwise
        labels = struct('da', 'd_a [nm]', 'dpp', 'd_pp [nm]', ...
            'dpp_individual', 'd_pp^(i,j) [nm]', ...
            'sigmapp', 'sigma_pp [-]', 'ca', 'c_a [-]', ...
            'zbar', 'z_a [-]', 'sbar', 's_a [-]', ...
            'n_subagg', 'N_subagg [-]', ...
            'f_col_subagg', 'f_col_subagg [-]');
end
label = labels.(quantity);

end

function state = logical_to_on_off(value)

if value
    state = 'on';
else
    state = 'off';
end

end

function rgb = color_list_to_rgb(colors)

colors = cellstr(colors);
rgb = zeros(numel(colors), 3);
for i = 1:numel(colors)
    rgb(i, :) = hex_to_rgb(colors{i});
end

end

function mask = interval_mask(values, interval)

if interval.include_lower
    lower_mask = values >= interval.lower;
else
    lower_mask = values > interval.lower;
end
if interval.include_upper
    upper_mask = values <= interval.upper;
else
    upper_mask = values < interval.upper;
end
mask = lower_mask & upper_mask;

end

function values = collapse_fraction_values(aggregate_table, entry_id)

% The collapse fraction is defined for each aggregate as the collapsed count
% divided by its subaggregate count. Invalid counts are rejected explicitly.

n_subagg = double(aggregate_table.n_subagg);
n_collapsed = double(aggregate_table.n_colaps);
mask = isfinite(n_subagg) & n_subagg >= 1 & isfinite(n_collapsed);
if any(n_collapsed(mask) < 0 | n_collapsed(mask) > n_subagg(mask))
    error('PFAL:main_tem_analysis_v1:InvalidCollapseCounts', ...
        ['Entry "%s" contains collapse counts outside the interval ', ...
        '[0, n_subagg].'], entry_id);
end
values = n_collapsed(mask) ./ n_subagg(mask);

end

function [matrix, valid_n] = frequency_matrix(frequency_table, ...
    entry_results, bins)

matrix = zeros(numel(entry_results), numel(bins));
valid_n = zeros(numel(entry_results), 1);
for i = 1:numel(entry_results)
    entry_id = entry_results(i).entry.id;
    for j = 1:numel(bins)
        row = strcmp(frequency_table.entry_id, entry_id) & ...
            strcmp(frequency_table.bin_id, bins(j).id);
        if nnz(row) ~= 1
            error('PFAL:main_tem_analysis_v1:MissingFrequencyRow', ...
                'Expected one frequency row for condition "%s", bin "%s".', ...
                entry_id, bins(j).id);
        end
        matrix(i, j) = frequency_table.frequency_percent(row);
        valid_n(i) = frequency_table.valid_n(row);
    end
end

end

function idx = weighted_sample_indices(n_values, n_samples, weights)

weights = double(weights(:));
weights(~isfinite(weights) | weights < 0) = 0;
if sum(weights) <= 0
    weights = ones(n_values, 1);
end
edges = cumsum(weights / sum(weights));
edges(end) = 1;
r = rand(n_samples, 1);
idx = arrayfun(@(x) find(edges >= x, 1, 'first'), r);

end

function [support, matrix] = exact_count_frequency_matrix( ...
    frequency_table, entry_results)

support = unique(frequency_table.n_subagg, 'sorted');
matrix = zeros(numel(support), numel(entry_results));
for i = 1:numel(entry_results)
    entry_id = entry_results(i).entry.id;
    rows = strcmp(frequency_table.entry_id, entry_id);
    entry_support = frequency_table.n_subagg(rows);
    entry_frequency = frequency_table.frequency_percent(rows);
    [entry_support, order] = sort(entry_support);
    entry_frequency = entry_frequency(order);
    if ~isequal(entry_support, support)
        error('PFAL:main_tem_analysis_v1:IncompleteExactFrequency', ...
            ['The exact-frequency summary for condition "%s" does not ', ...
            'cover the common integer support.'], entry_id);
    end
    matrix(:, i) = entry_frequency;
end

end

function lim = configured_axis_limit(cfg_tem, figure_name, field_name, ...
    default_lim, data_values, scale_type)

override = configured_plot_value(cfg_tem, figure_name, field_name);
if ~isempty(override)
    lim = double(override(:)).';
    if numel(lim) ~= 2 || any(~isfinite(lim)) || lim(2) <= lim(1)
        error('PFAL:main_tem_analysis_v1:InvalidAxisLimit', ...
            'plots.axis_overrides.%s.%s must be a finite [min max] vector.', ...
            figure_name, field_name);
    end
    return
end

lim = expand_limits_to_data(default_lim, data_values, scale_type);

end

function ticks = configured_axis_ticks(cfg_tem, figure_name, field_name, ...
    default_ticks, current_lim, scale_type)

override = configured_plot_value(cfg_tem, figure_name, field_name);
if ~isempty(override)
    ticks = double(override(:)).';
    return
end

ticks = default_ticks(default_ticks >= current_lim(1) & ...
    default_ticks <= current_lim(2));
if isempty(ticks) && strcmp(scale_type, 'log')
    ticks = log_ticks_for_limits(current_lim);
elseif isempty(ticks)
    ticks = linspace(current_lim(1), current_lim(2), 5);
end

end

function value = configured_plot_value(cfg_tem, figure_name, field_name)

value = [];
if ~isfield(cfg_tem.plots.figures, figure_name)
    return
end
figure_cfg = cfg_tem.plots.figures.(figure_name);
field_map = struct('xlim', 'x_limits', 'ylim', 'y_limits', ...
    'xticks', 'x_ticks', 'yticks', 'y_ticks');
if isfield(field_map, field_name)
    canonical_name = field_map.(field_name);
    if isfield(figure_cfg, canonical_name) && ...
            ~isempty(figure_cfg.(canonical_name))
        value = figure_cfg.(canonical_name);
        return
    end
end
if isfield(figure_cfg, 'axis_overrides') && ...
        isstruct(figure_cfg.axis_overrides) && ...
        isfield(figure_cfg.axis_overrides, field_name) && ...
        ~isempty(figure_cfg.axis_overrides.(field_name))
    value = figure_cfg.axis_overrides.(field_name);
end

end

function lim = expand_limits_to_data(default_lim, data_values, scale_type)

data_values = data_values(isfinite(data_values));
if strcmp(scale_type, 'log')
    data_values = data_values(data_values > 0);
end
if isempty(data_values)
    lim = default_lim;
    return
end

if strcmp(scale_type, 'log')
    if min(data_values) < default_lim(1) || max(data_values) > default_lim(2)
        lim = nice_log_limits(data_values, [0.9 1.1]);
        lim = [min(lim(1), default_lim(1)), max(lim(2), default_lim(2))];
    else
        lim = default_lim;
    end
else
    data_span = max(data_values) - min(data_values);
    pad = max(0.05 * data_span, 0.05 * max(abs(data_values)));
    if pad == 0
        pad = 0.05;
    end
    lim = [min(default_lim(1), min(data_values) - pad), ...
        max(default_lim(2), max(data_values) + pad)];
end

end

function lim = nice_log_limits(values, factors)

values = values(isfinite(values) & values > 0);
if isempty(values)
    lim = [1 10];
    return
end

raw_lim = [min(values) * factors(1), max(values) * factors(2)];
lim = [10 ^ floor(log10(raw_lim(1))), 10 ^ ceil(log10(raw_lim(2)))];
if lim(1) < raw_lim(1) / 3
    lim(1) = floor(raw_lim(1) / 10 ^ floor(log10(raw_lim(1)))) * ...
        10 ^ floor(log10(raw_lim(1)));
end
if lim(2) > raw_lim(2) * 3
    lim(2) = ceil(raw_lim(2) / 10 ^ floor(log10(raw_lim(2)))) * ...
        10 ^ floor(log10(raw_lim(2)));
end
lim(1) = max(lim(1), realmin);

end

function lim = nice_linear_limits(values, padding_fraction, bounds)

values = values(isfinite(values));
if isempty(values)
    if isempty(bounds)
        lim = [0 1];
    else
        lim = bounds;
    end
    return
end
span_value = max(values) - min(values);
pad = padding_fraction * span_value;
if pad == 0
    pad = max(0.05 * abs(max(values)), 0.05);
end
lim = [min(values) - pad, max(values) + pad];
if ~isempty(bounds)
    lim(1) = max(bounds(1), lim(1));
    lim(2) = min(bounds(2), lim(2));
end
if lim(2) <= lim(1)
    lim = lim + [-0.5 0.5] * max(abs(lim(1)), 1) * 0.1;
end

end

function ticks = log_ticks_for_limits(lim)

lo_exp = floor(log10(lim(1)));
hi_exp = ceil(log10(lim(2)));
ticks = [];
for exponent = lo_exp:hi_exp
    ticks = [ticks, [1 2 5] * 10 ^ exponent]; %#ok<AGROW>
end
ticks = ticks(ticks >= lim(1) & ticks <= lim(2));

end

function text_out = latex_text(text_in)

text_out = char(text_in);
text_out = strrep(text_out, '\', '\\');
text_out = strrep(text_out, '_', '\_');
text_out = strrep(text_out, '%', '\%');

end

function text_out = tex_plain_text(text_in)

% Escape plain legend text without inserting TeX font-family commands. This
% keeps the legend in the same configured typeface as axes and labels while
% still protecting underscores, percent signs, and literal backslashes.
text_out = latex_text(text_in);

end

function export_figure(fig, cfg_tem, name)

% Direct raster export remains available for minimal installations. The paper
% profile renders its PNG from the final vector PDF because that route preserves
% the embedded typeface, text spacing, colors, and geometry exactly.

if cfg_tem.plots.export
    pdf_path = fullfile(cfg_tem.outputs.results_root, sprintf('%s.pdf', name));
    png_path = fullfile(cfg_tem.outputs.results_root, sprintf('%s.png', name));
    if cfg_tem.plots.png && strcmp(cfg_tem.plots.png_source, 'figure')
        finalize_figure_text(fig, cfg_tem)
        exportgraphics(fig, png_path, 'BackgroundColor', 'white', ...
            'Resolution', cfg_tem.plots.png_resolution)
    end
    if cfg_tem.plots.pdf
        finalize_figure_text(fig, cfg_tem)
        exportgraphics(fig, pdf_path, 'ContentType', 'vector', ...
            'BackgroundColor', 'white')
    end
    if cfg_tem.plots.png && strcmp(cfg_tem.plots.png_source, 'pdf')
        UTILS.RASTERIZE_PDF(pdf_path, png_path, ...
            cfg_tem.plots.png_resolution)
    end
    if cfg_tem.plots.save_figures
        finalize_figure_text(fig, cfg_tem)
        savefig(fig, fullfile(cfg_tem.outputs.results_root, ...
            sprintf('%s.fig', name)))
    end
end
if strcmpi(cfg_tem.plots.visible, 'on')
    fig.Visible = 'on';
    drawnow
end

end

function finalize_figure_text(fig, cfg_tem)

% Reapply the configured typeface and weight after all plot and legend
% objects exist because several MATLAB chart constructors create child text
% objects with their own defaults. Normalized label offsets keep scientific
% axis labels at a consistent distance from every plot frame.
font_name = resolve_font(cfg_tem.plots.font);
axes_handles = findall(fig, 'Type', 'axes');
for i = 1:numel(axes_handles)
    ax = axes_handles(i);
    set(ax, 'FontName', font_name, ...
        'FontWeight', cfg_tem.plots.font.weight, ...
        'XTickLabelRotation', cfg_tem.plots.axes.tick_label_rotation)
    if strcmpi(cfg_tem.plots.axes.label_position_mode, 'manual')
        position_axis_label(ax.XLabel, [0.5, ...
            -cfg_tem.plots.axes.x_label_offset])
        position_axis_label(ax.YLabel, [...
            -cfg_tem.plots.axes.y_label_offset, 0.5])
    end
    set(ax.Title, 'FontName', font_name, ...
        'FontWeight', cfg_tem.plots.font.weight)
end

text_handles = findall(fig, 'Type', 'text');
set(text_handles, 'FontName', font_name, ...
    'FontWeight', cfg_tem.plots.font.weight)

legend_handles = findall(fig, 'Type', 'legend');
set(legend_handles, 'FontName', font_name, ...
    'FontWeight', cfg_tem.plots.font.weight, ...
    'FontAngle', 'normal')
drawnow

end

function position_axis_label(label_handle, normalized_xy)

% Empty labels retain MATLAB's default position so auxiliary axes and chart
% internals are not moved unnecessarily.
if isempty(label_handle.String)
    return
end
label_handle.Units = 'normalized';
position = label_handle.Position;
position(1:2) = normalized_xy;
label_handle.Position = position;

end

function write_optional_table(table_in, file_path)

if ~isempty(table_in)
    writetable(table_in, file_path)
end

end

function table_out = vertcat_nonempty(table_cells)

table_cells = table_cells(~cellfun(@isempty, table_cells));
if isempty(table_cells)
    table_out = table();
else
    table_out = vertcat(table_cells{:});
end

end

function value = require_agg_numeric_field(Aggs, idx, field_name)

if ~isfield(Aggs, field_name)
    error('PFAL:main_tem_analysis_v1:MissingAggregateField', ...
        'Aggs is missing required field "%s".', field_name);
end
value = double(Aggs(idx).(field_name));
if ~isscalar(value) || ~isfinite(value)
    error('PFAL:main_tem_analysis_v1:InvalidAggregateField', ...
        'Aggs(%d).%s must be a finite scalar.', idx, field_name);
end

end

function value = optional_agg_numeric_field(Aggs, idx, field_name)

value = NaN;
if isfield(Aggs, field_name) && ~isempty(Aggs(idx).(field_name))
    value = double(Aggs(idx).(field_name));
    if ~isscalar(value) || ~isfinite(value)
        value = NaN;
    end
end

end

function rgb = hex_to_rgb(hex_value)

hex_value = char(hex_value);
if startsWith(hex_value, '#')
    hex_value = hex_value(2:end);
end
if numel(hex_value) ~= 6
    error('PFAL:main_tem_analysis_v1:InvalidHexColor', ...
        'Hex color must contain six characters.');
end
rgb = [hex2dec(hex_value(1:2)), hex2dec(hex_value(3:4)), ...
    hex2dec(hex_value(5:6))] / 255;

end
