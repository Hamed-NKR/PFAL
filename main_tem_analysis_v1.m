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

%% Process configured TEM entries

n_entries = numel(cfg_tem.entries);
entry_result_cells = cell(n_entries, 1);
aggregate_tables = cell(n_entries, 1);
primary_tables = cell(n_entries, 1);
summary_tables = cell(n_entries, 1);
model_input_tables = cell(n_entries, 1);

for i = 1:n_entries
    entry_result_cells{i} = process_tem_entry(cfg_tem.entries(i), cfg_tem.analysis);
    aggregate_tables{i} = entry_result_cells{i}.aggregate_table;
    primary_tables{i} = entry_result_cells{i}.primary_table;
    summary_tables{i} = entry_result_cells{i}.summary_table;
    model_input_tables{i} = entry_result_cells{i}.model_inputs;
end

entry_results = vertcat(entry_result_cells{:});
aggregate_table = vertcat_nonempty(aggregate_tables);
primary_table = vertcat_nonempty(primary_tables);
summary_table = vertcat_nonempty(summary_tables);
model_inputs = vertcat_nonempty(model_input_tables);

% Build the figure-facing summary products before saving.  These tables are
% intentionally separate from model_inputs: model_inputs contains the four
% scalar quantities consumed by main_scale/main_scatter, while the plot
% artifacts preserve the publication figure quantities from the ATEMS
% workflow for audit, reuse, and visual regression checks.
plot_artifacts = build_plot_artifacts(entry_results, cfg_tem);
ensemble_primary_summary = plot_artifacts.ensemble_primary_summary;
per_aggregate_dpp_summary = plot_artifacts.per_aggregate_dpp_summary;
per_aggregate_sigmapp_summary = plot_artifacts.per_aggregate_sigmapp_summary;
hybridity_frequency_summary = plot_artifacts.hybridity_frequency_summary;
collapse_frequency_summary = plot_artifacts.collapse_frequency_summary;
plot_diagnostics = plot_artifacts.plot_diagnostics;

%% Export machine-readable outputs

metadata = struct( ...
    'generated_at', char(datetime('now')), ...
    'config_file', cfg_tem.config_file, ...
    'entry_count', n_entries);

save(fullfile(cfg_tem.outputs.data_root, cfg_tem.outputs.mat_file), ...
    'cfg_tem', 'metadata', 'entry_results', 'aggregate_table', ...
    'primary_table', 'summary_table', 'model_inputs', 'plot_artifacts', ...
    'ensemble_primary_summary', 'per_aggregate_dpp_summary', ...
    'per_aggregate_sigmapp_summary', 'hybridity_frequency_summary', ...
    'collapse_frequency_summary', 'plot_diagnostics', '-v7.3')

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
    aggregate_ids = (1:numel(Aggs)).';
else
    aggregate_ids = entry.aggregate_ids(:);
end

if any(aggregate_ids > numel(Aggs))
    error('PFAL:main_tem_analysis_v1:AggregateIdOutOfRange', ...
        'Entry "%s" includes aggregate ids larger than numel(Aggs) = %d.', ...
        entry.id, numel(Aggs));
end

aggregate_rows = table();
primary_rows = table();

for j = 1:numel(aggregate_ids)
    agg_id = aggregate_ids(j);
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
    % aggregate projected area.  ATEMS used this value to decide when manual
    % sizing likely under-sampled larger primary particles.  PFAL keeps the
    % same smooth logistic switch: low-coverage aggregates receive an
    % inverse-area correction, while high-coverage aggregates approach
    % unweighted statistics continuously instead of through a hard cutoff.
    coverage = sum(dpp_values .^ 2) / (da_value ^ 2);
    alpha_value = 1 - 1 / (1 + exp( ...
        -(coverage - analysis_cfg.coverage_threshold) / ...
        analysis_cfg.logistic_bandwidth));

    % Per-aggregate weights estimate the geometric mean and GSD of primary
    % particles within this one aggregate.  Ensemble weights below include
    % parent aggregate area as in ATEMS because the ensemble distribution is
    % interpreted as a condition-level primary-particle population, not just
    % a concatenation of equal-weight manual measurements.
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

summary_table = summarize_entry(entry, aggregate_rows, primary_rows, ...
    analysis_cfg);
model_inputs = build_model_inputs(entry, aggregate_rows);

entry_result = struct();
entry_result.entry = entry;
entry_result.aggregate_table = aggregate_rows;
entry_result.primary_table = primary_rows;
entry_result.summary_table = summary_table;
entry_result.model_inputs = model_inputs;

end

function summary_table = summarize_entry(entry, aggregate_table, ...
    primary_table, analysis_cfg)

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
    'VariableNames', {'entry_id', 'entry_label', 'n_aggregates', ...
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
    build_hybridity_frequency_summary(entry_results);
plot_artifacts.collapse_frequency_summary = ...
    build_collapse_frequency_summary(entry_results);
plot_artifacts.plot_diagnostics = build_plot_diagnostics( ...
    entry_results, entry_styles, cfg_tem);

end

function plot_tem_analysis(entry_results, cfg_tem, plot_artifacts)

% The ATEMS script created a coherent publication figure set rather than
% one-off diagnostic plots.  Keep that structure here and make the entry
% selection data-driven: each configured entry can participate in the
% publication plots, and each figure decides whether its required metrics are
% present before drawing.
plot_dpp_vs_da_publication(entry_results, cfg_tem, plot_artifacts)
plot_primary_particle_distributions(entry_results, cfg_tem, plot_artifacts)
plot_aggregate_metric_distributions(entry_results, cfg_tem, plot_artifacts)
plot_hybridity_collapse_frequencies(entry_results, cfg_tem, plot_artifacts)
plot_dpp_vs_da_by_hybridity(entry_results, cfg_tem, plot_artifacts)

end

function plot_dpp_vs_da_publication(entry_results, cfg_tem, plot_artifacts)

[entry_results, styles] = select_publication_entries(entry_results, ...
    plot_artifacts.entry_styles);
assert_required_plot_metric(entry_results, 'da_nm', 'dpp_vs_da_manual_tem')
assert_required_plot_metric(entry_results, 'dbarpp_nm', 'dpp_vs_da_manual_tem')

fig = figure('Visible', cfg_tem.plots.visible);
fig.Position = publication_position('dpp_vs_da', numel(entry_results));
set(fig, 'color', 'white')

% Reference curve from Olfert and Rogak (2019), plotted exactly as in the
% ATEMS workflow.  It anchors the manual TEM measurements to the universal
% projected-area-diameter/primary-diameter relation used elsewhere in PFAL.
[da0, dpp0] = reference_dpp_da_curve();
plt_ref = plot(da0, dpp0, 'Color', hex_to_rgb('#DEAA79'), ...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

entry_handles = gobjects(numel(entry_results), 1);
legend_entries = cell(numel(entry_results), 1);
all_da = [];
all_dpp = [];
for i = 1:numel(entry_results)
    agg_table = entry_results(i).aggregate_table;
    style = styles(i);
    marker_size = style.scatter_size;
    entry_handles(i) = scatter(agg_table.da_nm, agg_table.dbarpp_nm, ...
        marker_size, hex_to_rgb(style.scatter_color), style.marker, ...
        'LineWidth', 1.5);
    legend_entries{i} = sprintf('%s (n = %d)', ...
        latex_text(entry_results(i).entry.label), height(agg_table));
    all_da = [all_da; agg_table.da_nm(:)]; %#ok<AGROW>
    all_dpp = [all_dpp; agg_table.dbarpp_nm(:)]; %#ok<AGROW>
end

ax = gca;
set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', 16, ...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim(configured_axis_limit(cfg_tem, 'dpp_vs_da_manual_tem', 'xlim', ...
    nice_log_limits(all_da, [0.85 1.2]), all_da, 'log'))
ylim(configured_axis_limit(cfg_tem, 'dpp_vs_da_manual_tem', 'ylim', ...
    nice_log_limits(all_dpp, [0.95 1.05]), all_dpp, 'log'))
xticks(configured_axis_ticks(cfg_tem, 'dpp_vs_da_manual_tem', 'xticks', ...
    atems_da_ticks(), xlim(ax), 'log'))
xtickangle(90)
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 24)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 24)
box on
legend([entry_handles; plt_ref], ...
    [legend_entries; {'Olfert and Rogak (2019)'}], ...
    'interpreter', 'latex', 'FontSize', 16, ...
    'location', 'southoutside', 'Orientation', 'horizontal', ...
    'NumColumns', min(2, numel(entry_results) + 1))

export_figure(fig, cfg_tem, 'dpp_vs_da_manual_tem')

end

function plot_primary_particle_distributions(entry_results, cfg_tem, ...
    plot_artifacts)

[entry_results, styles] = select_publication_entries(entry_results, ...
    plot_artifacts.entry_styles);
assert_required_plot_metric(entry_results, 'dpp_nm', ...
    'primary_particle_distributions')
assert_required_plot_metric(entry_results, 'dbarpp_nm', ...
    'primary_particle_distributions')
assert_required_plot_metric(entry_results, 'sigmapp', ...
    'primary_particle_distributions')

fig = figure('Visible', cfg_tem.plots.visible);
fig.Position = publication_position('primary_particle_distributions', ...
    numel(entry_results));
set(fig, 'color', 'white')
tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact')

% Panel 1 follows the ATEMS weighted-ensemble logic.  Individual ImageJ
% primary-particle diameters are weighted by the aggregate-to-primary area
% ratio with a logistic coverage switch.  The resampled boxplot makes the
% weighted ensemble visible while the KDE ridge preserves the continuous
% distribution shape.
nexttile
plot_weighted_primary_box_and_kde(entry_results, styles, cfg_tem)
all_dpp = collect_primary_metric(entry_results, 'dpp_nm');
ylim(configured_axis_limit(cfg_tem, 'primary_particle_distributions', ...
    'ensemble_dpp_ylim', [5 55], all_dpp, 'log'))
yticks(configured_axis_ticks(cfg_tem, 'primary_particle_distributions', ...
    'ensemble_dpp_yticks', [5 10 20 40 80], ylim, 'log'))

% Panel 2 compares aggregate-level geometric means of manually sized
% primaries.  This is the quantity later summarized as GM_dpp/GSD_dpp for
% PFAL model inputs, so it stays visually distinct from the raw ensemble.
nexttile
plot_publication_box_metric(entry_results, styles, 'dbarpp_nm', ...
    '$d_\mathrm{pp}$ [nm]', 0.3)
all_dbarpp = collect_aggregate_metric(entry_results, 'dbarpp_nm');
ylim(configured_axis_limit(cfg_tem, 'primary_particle_distributions', ...
    'aggregate_mean_dpp_ylim', [8 24], all_dbarpp, 'linear'))

% Panel 3 compares the within-aggregate geometric standard deviation of
% primary-particle diameter.  It is not a model input itself, but it is the
% primary check that the manual sizing spread was carried over correctly.
nexttile
plot_publication_box_metric(entry_results, styles, 'sigmapp', ...
    '$\sigma_\mathrm{pp}$ [-]', 0.3)
all_sigmapp = collect_aggregate_metric(entry_results, 'sigmapp');
ylim(configured_axis_limit(cfg_tem, 'primary_particle_distributions', ...
    'sigmapp_ylim', [1.11 1.62], all_sigmapp, 'linear'))
yticks(configured_axis_ticks(cfg_tem, 'primary_particle_distributions', ...
    'sigmapp_yticks', [1.2 1.3 1.4 1.5 1.6], ylim, 'linear'))

export_figure(fig, cfg_tem, 'primary_particle_distributions')

end

function plot_aggregate_metric_distributions(entry_results, cfg_tem, ...
    plot_artifacts)

[entry_results, styles] = select_publication_entries(entry_results, ...
    plot_artifacts.entry_styles);
assert_required_plot_metric(entry_results, 'da_nm', ...
    'aggregate_metric_distributions')

fig = figure('Visible', cfg_tem.plots.visible);
fig.Position = publication_position('aggregate_metric_distributions', ...
    numel(entry_results));
set(fig, 'color', 'white')
tiledlayout(2, 2, 'Padding', 'loose', 'TileSpacing', 'compact')

% ATEMS used this 2x2 figure to separate aggregate size from morphology.
% When n_subagg is available, rows with n_subagg <= 0 are excluded from
% morphology distributions because those aggregates have no meaningful
% hybrid-subaggregate segmentation.
nexttile
plot_publication_box_metric(entry_results, styles, 'da_nm', ...
    '$d_\mathrm{a}$ [nm]', 0.25, true)
set(gca, 'YScale', 'log')
all_da = collect_aggregate_metric(entry_results, 'da_nm', true);
ylim(configured_axis_limit(cfg_tem, 'aggregate_metric_distributions', ...
    'da_ylim', nice_log_limits(all_da, [0.85 1.15]), all_da, 'log'))

nexttile
plot_publication_box_metric(entry_results, styles, 'ca', ...
    '$c_\mathrm{a}$ [-]', 0.25, true)
all_ca = collect_aggregate_metric(entry_results, 'ca', true);
ylim(configured_axis_limit(cfg_tem, 'aggregate_metric_distributions', ...
    'ca_ylim', nice_linear_limits(all_ca, 0.08, [0 1]), all_ca, 'linear'))

nexttile
plot_publication_box_metric(entry_results, styles, 'zbar_opt', ...
    '$z_\mathrm{a}$ [-]', 0.25, true)
all_zbar = collect_aggregate_metric(entry_results, 'zbar_opt', true);
ylim(configured_axis_limit(cfg_tem, 'aggregate_metric_distributions', ...
    'zbar_opt_ylim', nice_linear_limits(all_zbar, 0.08, []), all_zbar, ...
    'linear'))

nexttile
plot_publication_box_metric(entry_results, styles, 'sbar_opt', ...
    '$s_\mathrm{a}$ [-]', 0.25, true)
all_sbar = collect_aggregate_metric(entry_results, 'sbar_opt', true);
ylim(configured_axis_limit(cfg_tem, 'aggregate_metric_distributions', ...
    'sbar_opt_ylim', nice_linear_limits(all_sbar, 0.08, []), all_sbar, ...
    'linear'))

export_figure(fig, cfg_tem, 'aggregate_metric_distributions')

end

function plot_hybridity_collapse_frequencies(entry_results, cfg_tem, ...
    plot_artifacts)

[entry_results, ~, mask] = select_publication_entries(entry_results, ...
    plot_artifacts.entry_styles);
assert_required_plot_metric(entry_results, 'n_subagg', ...
    'hybridity_collapse_frequencies')
assert_required_plot_metric(entry_results, 'n_colaps', ...
    'hybridity_collapse_frequencies')

hyb_table = plot_artifacts.hybridity_frequency_summary(mask, :);
clp_table = plot_artifacts.collapse_frequency_summary(mask, :);

fig = figure('Visible', cfg_tem.plots.visible);
fig.Position = publication_position('hybridity_collapse_frequencies', ...
    numel(entry_results));
set(fig, 'color', 'white')
tiledlayout(1, 2, 'Padding', 'loose', 'TileSpacing', 'compact')

nexttile
hyb_matrix = [ ...
    hyb_table.frequency_percent_n_hyb_eq_1, ...
    hyb_table.frequency_percent_n_hyb_eq_2, ...
    hyb_table.frequency_percent_n_hyb_3_to_5, ...
    hyb_table.frequency_percent_n_hyb_6_to_10, ...
    hyb_table.frequency_percent_n_hyb_gt_10];
bar_handles = bar(hyb_matrix, 'stacked');
style_stacked_bars(bar_handles, atems_hybridity_colors())
format_frequency_axis(hyb_table.entry_label, hyb_table.valid_n)
ylabel('Frequency [$\%$]', 'interpreter', 'latex', 'FontSize', 22)
lgd = legend({'$n_\mathrm{hyb} = 1$', '$n_\mathrm{hyb} = 2$', ...
    '$3 \le n_\mathrm{hyb} \le 5$', ...
    '$6 \le n_\mathrm{hyb} \le 10$', '$n_\mathrm{hyb} > 10$'}, ...
    'interpreter', 'latex', 'FontSize', 18, 'location', ...
    'northoutside', 'orientation', 'horizontal', 'NumColumns', 2);
lgd.ItemTokenSize = [15, 15];
add_hybridity_arrows(gca, bar_handles, min(2, height(hyb_table)))

nexttile
collapse_matrix = [ ...
    clp_table.frequency_percent_r_clp_eq_0, ...
    clp_table.frequency_percent_r_clp_0_to_033, ...
    clp_table.frequency_percent_r_clp_033_to_067, ...
    clp_table.frequency_percent_r_clp_067_to_1, ...
    clp_table.frequency_percent_r_clp_eq_1];
bar_handles = bar(collapse_matrix, 'stacked');
style_stacked_bars(bar_handles, atems_collapse_colors())
format_frequency_axis(clp_table.entry_label, clp_table.valid_n)
lgd = legend({'$r_\mathrm{clp} = 0$', ...
    '$0 < r_\mathrm{clp} < 0.33$', ...
    '$0.33 \le r_\mathrm{clp} \le 0.67$', ...
    '$0.67 < r_\mathrm{clp} < 1$', '$r_\mathrm{clp} = 1$'}, ...
    'interpreter', 'latex', 'FontSize', 18, 'location', ...
    'northoutside', 'orientation', 'horizontal', 'NumColumns', 2);
lgd.ItemTokenSize = [15, 15];

export_figure(fig, cfg_tem, 'hybridity_collapse_frequencies')

end

function plot_dpp_vs_da_by_hybridity(entry_results, cfg_tem, plot_artifacts)

[entry_results, ~] = select_hybridity_scatter_entries(entry_results, ...
    plot_artifacts.entry_styles);
assert_required_plot_metric(entry_results, 'n_subagg', ...
    'dpp_vs_da_by_hybridity')
assert_required_plot_metric(entry_results, 'da_nm', ...
    'dpp_vs_da_by_hybridity')
assert_required_plot_metric(entry_results, 'dbarpp_nm', ...
    'dpp_vs_da_by_hybridity')

fig = figure('Visible', cfg_tem.plots.visible);
fig.Position = publication_position('dpp_vs_da_by_hybridity', ...
    numel(entry_results));
set(fig, 'color', 'white')

[da0, dpp0] = reference_dpp_da_curve();
plt_ref = plot(da0, dpp0, 'Color', hex_to_rgb('#DEAA79'), ...
    'LineStyle', '-.', 'LineWidth', 2.5);
hold on

all_da = collect_aggregate_metric(entry_results, 'da_nm');
all_dpp = collect_aggregate_metric(entry_results, 'dbarpp_nm');
all_n_hyb = collect_aggregate_metric(entry_results, 'n_subagg');

hyb_colors = atems_hybridity_colors();
plt_1 = scatter(all_da(all_n_hyb == 1), all_dpp(all_n_hyb == 1), ...
    25, hyb_colors(1, :), '^', 'LineWidth', 1.5);
plt_2 = scatter(all_da(all_n_hyb == 2), all_dpp(all_n_hyb == 2), ...
    35, hyb_colors(2, :), 's', 'LineWidth', 1.5);
mask_3_to_5 = all_n_hyb >= 3 & all_n_hyb <= 5;
plt_3 = scatter(all_da(mask_3_to_5), all_dpp(mask_3_to_5), ...
    35, hyb_colors(3, :), 'h', 'LineWidth', 1.5);
plt_gt_5 = scatter(all_da(all_n_hyb > 5), all_dpp(all_n_hyb > 5), ...
    25, hyb_colors(4, :), 'o', 'LineWidth', 1.5);

ax = gca;
box on
set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', 11, ...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim(configured_axis_limit(cfg_tem, 'dpp_vs_da_by_hybridity', 'xlim', ...
    nice_log_limits(all_da, [0.8 1.2]), all_da, 'log'))
ylim(configured_axis_limit(cfg_tem, 'dpp_vs_da_by_hybridity', 'ylim', ...
    nice_log_limits(all_dpp, [0.95 1.05]), all_dpp, 'log'))
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$\overline{d}_\mathrm{pp}$ [nm]', 'interpreter', 'latex', ...
    'FontSize', 14)
legend([plt_ref, plt_1, plt_2, plt_3, plt_gt_5], ...
    {'Olfert and Rogak (2019)', '$n_\mathrm{hyb} = 1$', ...
    '$n_\mathrm{hyb} = 2$', '$3 \le n_\mathrm{hyb} \le 5$', ...
    '$n_\mathrm{hyb} > 5$'}, 'interpreter', 'latex', ...
    'FontSize', 12, 'location', 'northoutside', ...
    'orientation', 'horizontal', 'NumColumns', 2)

export_figure(fig, cfg_tem, 'dpp_vs_da_by_hybridity')

end

function plot_weighted_primary_box_and_kde(entry_results, styles, cfg_tem)

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

boxplot(values, groups, 'Labels', labels, 'Notch', 'on', ...
    'Symbol', 'o', 'Widths', 0.25, 'Colors', box_line_colors(styles))
style_boxplot(gca, styles)
hold on

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
    x_anchor = i - 0.3;
    density_x = x_anchor - 0.25 * f / max(f);
    plot(density_x, y, 'Color', hex_to_rgb(styles(i).box_edge_color), ...
        'LineWidth', 1.25)
    fill([density_x, x_anchor * ones(size(density_x))], ...
        [y, fliplr(y)], hex_to_rgb(styles(i).box_face_color), ...
        'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

set(gca, 'TickLabelInterpreter', 'tex', 'FontSize', 16, ...
    'TickLength', [0.015 0.015], 'YScale', 'log')
ylabel('$d_\mathrm{pp}^\mathrm{(i)}$ [nm]', 'interpreter', ...
    'latex', 'FontSize', 24)
xlim([0.3, numel(entry_results) + 0.3])
box on

end

function plot_publication_box_metric(entry_results, styles, metric_name, ...
    y_label, box_width, varargin)

filter_hybrid_rows = false;
if nargin > 5
    filter_hybrid_rows = varargin{1};
end

values = [];
groups = [];
labels = condition_labels(entry_results, 'aggregate');
for i = 1:numel(entry_results)
    metric_values = entry_results(i).aggregate_table.(metric_name);
    if filter_hybrid_rows && ismember('n_subagg', ...
            entry_results(i).aggregate_table.Properties.VariableNames)
        hybrid_mask = entry_results(i).aggregate_table.n_subagg > 0 | ...
            ~isfinite(entry_results(i).aggregate_table.n_subagg);
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

boxplot(values, groups, 'Labels', labels, 'Notch', 'on', ...
    'Symbol', 'o', 'Widths', box_width, 'Colors', box_line_colors(styles))
style_boxplot(gca, styles)
set(gca, 'TickLabelInterpreter', 'tex', 'FontSize', 14, ...
    'TickLength', [0.02 0.02])
ylabel(y_label, 'interpreter', 'latex', 'FontSize', 18)
box on

end

function [entry_results_out, styles_out, mask] = select_publication_entries( ...
    entry_results, styles)

mask = true(numel(entry_results), 1);
for i = 1:numel(entry_results)
    if isfield(entry_results(i).entry, 'include_in_publication_plots')
        mask(i) = entry_results(i).entry.include_in_publication_plots;
    end
end
entry_results_out = entry_results(mask);
styles_out = styles(mask);
if isempty(entry_results_out)
    error('PFAL:main_tem_analysis_v1:NoPublicationEntries', ...
        ['No TEM entries are enabled for publication plotting. Set ', ...
        'include_in_publication_plots=true for at least one entry.'])
end

end

function [entry_results_out, styles_out, mask] = select_hybridity_scatter_entries( ...
    entry_results, styles)

mask = true(numel(entry_results), 1);
for i = 1:numel(entry_results)
    if isfield(entry_results(i).entry, 'include_in_hybridity_scatter')
        mask(i) = entry_results(i).entry.include_in_hybridity_scatter;
    end
end
entry_results_out = entry_results(mask);
styles_out = styles(mask);
if isempty(entry_results_out)
    error('PFAL:main_tem_analysis_v1:NoHybridityScatterEntries', ...
        ['No TEM entries are enabled for the dpp-vs-da hybridity ', ...
        'scatter plot. Set include_in_hybridity_scatter=true for at ', ...
        'least one entry.'])
end

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

function styles = resolve_entry_styles(entry_results)

styles = repmat(default_style(), numel(entry_results), 1);
palette_index = 0;
for i = 1:numel(entry_results)
    entry = entry_results(i).entry;
    semantic_style = semantic_entry_style(entry.entry_type);
    explicit_color = '';
    explicit_marker = '';
    if isfield(entry, 'plot')
        if isfield(entry.plot, 'color')
            explicit_color = char(entry.plot.color);
        end
        if isfield(entry.plot, 'marker')
            explicit_marker = char(entry.plot.marker);
        end
    end

    if ~isempty(semantic_style.scatter_color)
        style = semantic_style;
    else
        palette_index = palette_index + 1;
        style = extended_entry_style(palette_index);
    end

    if ~isempty(explicit_color)
        style.scatter_color = explicit_color;
        if isempty(semantic_style.scatter_color)
            scatter_rgb = hex_to_rgb(explicit_color);
            style.box_face_color = rgb_to_hex(mix_rgb(scatter_rgb, ...
                [1 1 1], 0.35));
            style.box_edge_color = rgb_to_hex(mix_rgb(scatter_rgb, ...
                [0 0 0], 0.20));
            style.median_color = rgb_to_hex(mix_rgb(scatter_rgb, ...
                [0 0 0], 0.55));
        end
    end
    if ~isempty(explicit_marker)
        style.marker = explicit_marker;
    end

    style.entry_id = entry.id;
    style.entry_label = entry.label;
    style.entry_type = entry.entry_type;
    styles(i) = style;
end

end

function style = semantic_entry_style(entry_type)

style = default_style();
switch lower(char(entry_type))
    case {'low_agglomeration', 'low_agglom', 'lal'}
        style.scatter_color = '#C96868';
        style.box_face_color = '#DC8686';
        style.box_edge_color = '#8D493A';
        style.median_color = '#632626';
        style.marker = '^';
        style.scatter_size = 25;
    case {'moderate_collapse', 'moderate_agglomeration', ...
            'mod_collapse', 'hal'}
        style.scatter_color = '#006989';
        style.box_face_color = '#7EACB5';
        style.box_edge_color = '#537188';
        style.median_color = '#374259';
        style.marker = 'o';
        style.scatter_size = 25;
    case {'high_agglomeration', 'high_agglom', 'ex_agglomeration', ...
            'exaglom'}
        style.scatter_color = '#8174A0';
        style.box_face_color = '#A888B5';
        style.box_edge_color = '#8174A0';
        style.median_color = '#574964';
        style.marker = 's';
        style.scatter_size = 35;
    case {'extra_collapse', 'extreme_collapse', 'ex_collapse', ...
            'excolaps'}
        style.scatter_color = '#659287';
        style.box_face_color = '#B1C29E';
        style.box_edge_color = '#659287';
        style.median_color = '#5F6F65';
        style.marker = 'd';
        style.scatter_size = 35;
end

end

function style = extended_entry_style(index_value)

seed = {
    '#C96868', '#DC8686', '#8D493A', '#632626', '^', 25;
    '#006989', '#7EACB5', '#537188', '#374259', 'o', 25;
    '#8174A0', '#A888B5', '#8174A0', '#574964', 's', 35;
    '#659287', '#B1C29E', '#659287', '#5F6F65', 'd', 35};
markers = {'^', 'o', 's', 'd', 'h', 'v', 'p', '>'};
seed_index = mod(index_value - 1, size(seed, 1)) + 1;
cycle = floor((index_value - 1) / size(seed, 1));

style = default_style();
if cycle == 0
    style.scatter_color = seed{seed_index, 1};
    style.box_face_color = seed{seed_index, 2};
    style.box_edge_color = seed{seed_index, 3};
    style.median_color = seed{seed_index, 4};
else
    base_rgb = hex_to_rgb(seed{seed_index, 1});
    shift = min(0.18 + 0.08 * floor((cycle - 1) / 2), 0.45);
    if mod(cycle, 2) == 1
        scatter_rgb = mix_rgb(base_rgb, [1 1 1], shift);
    else
        scatter_rgb = mix_rgb(base_rgb, [0 0 0], shift);
    end
    style.scatter_color = rgb_to_hex(scatter_rgb);
    style.box_face_color = rgb_to_hex(mix_rgb(scatter_rgb, [1 1 1], 0.35));
    style.box_edge_color = rgb_to_hex(mix_rgb(scatter_rgb, [0 0 0], 0.20));
    style.median_color = rgb_to_hex(mix_rgb(scatter_rgb, [0 0 0], 0.55));
end
style.marker = markers{mod(index_value - 1, numel(markers)) + 1};
style.scatter_size = seed{seed_index, 6};

end

function style = default_style()

style = struct( ...
    'entry_id', '', ...
    'entry_label', '', ...
    'entry_type', '', ...
    'scatter_color', '', ...
    'box_face_color', '', ...
    'box_edge_color', '', ...
    'median_color', '', ...
    'marker', 'o', ...
    'scatter_size', 30);

end

function style_boxplot(ax, styles)

boxes = findobj(ax, 'Tag', 'Box');
for i = 1:numel(boxes)
    group_index = nearest_group_index(boxes(i), numel(styles));
    patch(get(boxes(i), 'XData'), get(boxes(i), 'YData'), ...
        hex_to_rgb(styles(group_index).box_face_color), ...
        'EdgeColor', hex_to_rgb(styles(group_index).box_edge_color), ...
        'FaceAlpha', 0.5);
end

medians = findobj(ax, 'Tag', 'Median');
for i = 1:numel(medians)
    group_index = nearest_group_index(medians(i), numel(styles));
    set(medians(i), 'Color', hex_to_rgb(styles(group_index).median_color), ...
        'LineWidth', 2);
end

outliers = findobj(ax, 'Tag', 'Outliers');
for i = 1:numel(outliers)
    group_index = nearest_group_index(outliers(i), numel(styles));
    edge_color = hex_to_rgb(styles(group_index).box_edge_color);
    outliers(i).Color = edge_color;
    outliers(i).MarkerEdgeColor = edge_color;
    outliers(i).MarkerFaceColor = edge_color;
    outliers(i).MarkerSize = 3;
end

set(findobj(ax, 'type', 'line', 'tag', 'Upper Whisker'), ...
    'linestyle', '-')
set(findobj(ax, 'type', 'line', 'tag', 'Lower Whisker'), ...
    'linestyle', '-')

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

function frequency_table = build_hybridity_frequency_summary(entry_results)

entry_id = {};
entry_label = {};
valid_n = [];
counts = zeros(numel(entry_results), 5);
freq = zeros(numel(entry_results), 5);

for i = 1:numel(entry_results)
    values = entry_results(i).aggregate_table.n_subagg;
    values = values(isfinite(values) & values >= 1);
    entry_id{end + 1, 1} = entry_results(i).entry.id; %#ok<AGROW>
    entry_label{end + 1, 1} = entry_results(i).entry.label; %#ok<AGROW>
    valid_n(end + 1, 1) = numel(values); %#ok<AGROW>
    counts(i, :) = [ ...
        nnz(values == 1), ...
        nnz(values == 2), ...
        nnz(values >= 3 & values <= 5), ...
        nnz(values >= 6 & values <= 10), ...
        nnz(values > 10)];
    if valid_n(i) > 0
        freq(i, :) = 100 * counts(i, :) / valid_n(i);
    end
end

frequency_table = table(entry_id, entry_label, valid_n, ...
    counts(:, 1), counts(:, 2), counts(:, 3), counts(:, 4), counts(:, 5), ...
    freq(:, 1), freq(:, 2), freq(:, 3), freq(:, 4), freq(:, 5), ...
    'VariableNames', {'entry_id', 'entry_label', 'valid_n', ...
    'count_n_hyb_eq_1', 'count_n_hyb_eq_2', 'count_n_hyb_3_to_5', ...
    'count_n_hyb_6_to_10', 'count_n_hyb_gt_10', ...
    'frequency_percent_n_hyb_eq_1', ...
    'frequency_percent_n_hyb_eq_2', ...
    'frequency_percent_n_hyb_3_to_5', ...
    'frequency_percent_n_hyb_6_to_10', ...
    'frequency_percent_n_hyb_gt_10'});

end

function frequency_table = build_collapse_frequency_summary(entry_results)

entry_id = {};
entry_label = {};
valid_n = [];
counts = zeros(numel(entry_results), 5);
freq = zeros(numel(entry_results), 5);

for i = 1:numel(entry_results)
    n_subagg = entry_results(i).aggregate_table.n_subagg;
    n_colaps = entry_results(i).aggregate_table.n_colaps;
    mask = isfinite(n_subagg) & n_subagg >= 1 & isfinite(n_colaps);
    ratio = double(n_colaps(mask)) ./ double(n_subagg(mask));
    entry_id{end + 1, 1} = entry_results(i).entry.id; %#ok<AGROW>
    entry_label{end + 1, 1} = entry_results(i).entry.label; %#ok<AGROW>
    valid_n(end + 1, 1) = numel(ratio); %#ok<AGROW>
    counts(i, :) = [ ...
        nnz(ratio == 0), ...
        nnz(ratio > 0 & ratio < 0.33), ...
        nnz(ratio >= 0.33 & ratio <= 0.67), ...
        nnz(ratio > 0.67 & ratio < 1), ...
        nnz(ratio == 1)];
    if valid_n(i) > 0
        freq(i, :) = 100 * counts(i, :) / valid_n(i);
    end
end

frequency_table = table(entry_id, entry_label, valid_n, ...
    counts(:, 1), counts(:, 2), counts(:, 3), counts(:, 4), counts(:, 5), ...
    freq(:, 1), freq(:, 2), freq(:, 3), freq(:, 4), freq(:, 5), ...
    'VariableNames', {'entry_id', 'entry_label', 'valid_n', ...
    'count_r_clp_eq_0', 'count_r_clp_0_to_033', ...
    'count_r_clp_033_to_067', 'count_r_clp_067_to_1', ...
    'count_r_clp_eq_1', 'frequency_percent_r_clp_eq_0', ...
    'frequency_percent_r_clp_0_to_033', ...
    'frequency_percent_r_clp_033_to_067', ...
    'frequency_percent_r_clp_067_to_1', ...
    'frequency_percent_r_clp_eq_1'});

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
include_publication = [];
include_hybridity = [];
missing_optional_metrics = {};
diagnostic = {};
optional_metrics = {'n_subagg', 'n_colaps', 'ca', 'zbar_opt', 'sbar_opt'};

for i = 1:numel(entry_results)
    agg_table = entry_results(i).aggregate_table;
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
    include_publication(end + 1, 1) = ...
        entry_results(i).entry.include_in_publication_plots; %#ok<AGROW>
    include_hybridity(end + 1, 1) = ...
        entry_results(i).entry.include_in_hybridity_scatter; %#ok<AGROW>
    missing_optional_metrics{end + 1, 1} = strjoin(missing, ','); %#ok<AGROW>
    diagnostic{end + 1, 1} = sprintf( ...
        'ATEMS-derived style resolved; config=%s; resolution=%d', ...
        cfg_tem.config_file, cfg_tem.plots.resolution); %#ok<AGROW>
end

diagnostics = table(entry_id, entry_label, entry_type, scatter_color, ...
    box_face_color, box_edge_color, median_color, marker, ...
    include_publication, include_hybridity, missing_optional_metrics, ...
    diagnostic, 'VariableNames', {'entry_id', 'entry_label', ...
    'entry_type', 'scatter_color', 'box_face_color', 'box_edge_color', ...
    'median_color', 'marker', 'include_in_publication_plots', ...
    'include_in_hybridity_scatter', 'missing_optional_metrics', ...
    'diagnostic'});

end

function labels = condition_labels(entry_results, count_kind)

labels = cell(numel(entry_results), 1);
for i = 1:numel(entry_results)
    switch count_kind
        case 'primary'
            n_value = height(entry_results(i).primary_table);
        otherwise
            n_value = height(entry_results(i).aggregate_table);
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

% ATEMS used every 10 nm from 30-100 nm.  That is retained as the source
% pattern, but thinned here to prevent overlapping vertical labels when the
% number of configured entries changes the figure width.
ticks = [30 50 70 100 200 300 500 700];

end

function colors = atems_hybridity_colors()

colors = [
    hex_to_rgb('#C96868');
    hex_to_rgb('#FADFA1');
    hex_to_rgb('#7EACB5');
    hex_to_rgb('#8174A0');
    hex_to_rgb('#574964')];

end

function colors = atems_collapse_colors()

colors = [
    hex_to_rgb('#DEAA79');
    hex_to_rgb('#FFE6A9');
    hex_to_rgb('#B1C29E');
    hex_to_rgb('#659287');
    hex_to_rgb('#3A4D39')];

end

function style_stacked_bars(bar_handles, colors)

for i = 1:numel(bar_handles)
    bar_handles(i).BarWidth = 0.4;
    bar_handles(i).FaceColor = colors(i, :);
end

end

function format_frequency_axis(entry_labels, valid_n)

labels = cell(numel(entry_labels), 1);
for i = 1:numel(entry_labels)
    labels{i} = sprintf('%s\\newline(n = %d)', latex_text(entry_labels{i}), ...
        valid_n(i));
end
set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, ...
    'TickLabelInterpreter', 'tex', 'FontSize', frequency_label_font(labels), ...
    'TickLength', [0.02 0.02])
xtickangle(0)
yticks([0 20 40 60 80 100])
ylim([0 100])
xlim([0.5, numel(labels) + 0.5])
box on

end

function font_size = frequency_label_font(labels)

max_label_len = max(cellfun(@numel, labels));
if max_label_len > 26
    font_size = 14;
else
    font_size = 18;
end

end

function add_hybridity_arrows(ax, bar_handles, n_arrows)

if isempty(bar_handles) || n_arrows == 0
    return
end
for i = 1:n_arrows
    y_top = bar_handles(1).YData(i);
    if ~isfinite(y_top) || y_top <= 0
        continue
    end
    head_y = max(2, y_top - max(6, 0.2 * y_top));
    tail_y = min(98, head_y + 14);
    add_down_arrow(ax, i, tail_y, head_y)
end

end

function add_down_arrow(ax, x_value, y_start, y_end)

fig = ancestor(ax, 'figure');
old_units = ax.Units;
ax.Units = 'normalized';
ax_pos = ax.Position;
ax.Units = old_units;
xl = xlim(ax);
yl = ylim(ax);
x_norm = ax_pos(1) + (x_value - xl(1)) / diff(xl) * ax_pos(3);
y_start_norm = ax_pos(2) + (y_start - yl(1)) / diff(yl) * ax_pos(4);
y_end_norm = ax_pos(2) + (y_end - yl(1)) / diff(yl) * ax_pos(4);
annotation(fig, 'textarrow', [x_norm x_norm], ...
    [y_start_norm y_end_norm], 'String', '', 'HeadStyle', ...
    'plain', 'LineWidth', 2.0, 'Color', [0 0 0]);

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

function pos = publication_position(figure_name, n_entries)

switch figure_name
    case 'dpp_vs_da'
        pos = [50, 50, max(650, 140 * n_entries + 360), 800];
    case 'primary_particle_distributions'
        pos = [100, 0, max(850, 175 * n_entries + 330), 1300];
    case 'aggregate_metric_distributions'
        pos = [150, 50, max(900, 180 * n_entries + 260), 900];
    case 'hybridity_collapse_frequencies'
        pos = [200, 150, max(1200, 260 * n_entries + 420), 600];
    case 'dpp_vs_da_by_hybridity'
        pos = [250, 200, max(550, 90 * n_entries + 370), 600];
    otherwise
        pos = [100, 100, max(650, 200 * n_entries), 600];
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
if ~isfield(cfg_tem.plots, 'axis_overrides') || ...
        ~isstruct(cfg_tem.plots.axis_overrides)
    return
end
if ~isfield(cfg_tem.plots.axis_overrides, figure_name)
    return
end
figure_cfg = cfg_tem.plots.axis_overrides.(figure_name);
if isfield(figure_cfg, field_name) && ~isempty(figure_cfg.(field_name))
    value = figure_cfg.(field_name);
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

function rgb = mix_rgb(rgb_a, rgb_b, amount_b)

rgb = (1 - amount_b) * rgb_a + amount_b * rgb_b;
rgb = min(max(rgb, 0), 1);

end

function hex_value = rgb_to_hex(rgb)

rgb = round(255 * min(max(rgb, 0), 1));
hex_value = sprintf('#%02X%02X%02X', rgb(1), rgb(2), rgb(3));

end

function export_figure(fig, cfg_tem, name)

if cfg_tem.plots.export
    exportgraphics(fig, fullfile(cfg_tem.outputs.results_root, ...
        sprintf('%s.%s', name, cfg_tem.plots.format)), ...
        'BackgroundColor', 'white', 'Resolution', cfg_tem.plots.resolution)
    if cfg_tem.plots.save_figures
        savefig(fig, fullfile(cfg_tem.outputs.results_root, ...
            sprintf('%s.fig', name)))
    end
end

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
