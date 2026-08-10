function result = RUN_MAIN_VALIDATION(cfg)
%RUN_MAIN_VALIDATION Compare configured LD2 results with experimental data.
%   RESULT = UTILS.RUN_MAIN_VALIDATION(CFG) builds the two manuscript
%   validation figures, evaluates experimental observations against the
%   configured simulation fits, and writes reproducible run artifacts.

%% Load only the validation-ready variables

[ld2_source_data, source_info.ld2] = ...
    UTILS.LOAD_VALIDATION_SOURCE(cfg.sources.ld2);
[tem_source_data, source_info.tem] = ...
    UTILS.LOAD_VALIDATION_SOURCE(cfg.sources.tem);
[density_source_data, source_info.effective_density] = ...
    UTILS.LOAD_VALIDATION_SOURCE(cfg.sources.effective_density);

validate_source_payloads(ld2_source_data, tem_source_data, density_source_data);
font_name = resolve_font(cfg.figures.font);

%% Map named conditions to LD2 snapshots and experimental rows

conditions = prepare_conditions(cfg.conditions, ld2_source_data, ...
    tem_source_data.aggregate_table, ...
    density_source_data.effective_density_data, ...
    cfg.physics.material_density_kg_m3);
if isempty(conditions)
    error('PFAL:RUN_MAIN_VALIDATION:NoEnabledConditions', ...
        'No enabled validation conditions were found in the config.');
end

%% Build both validation panels and quantitative comparisons

[figure_dpp, dpp_metrics, dpp_predictions, dpp_fits] = ...
    plot_dpp_vs_da(conditions, cfg, font_name);
[figure_density, density_metrics, density_predictions, density_fits] = ...
    plot_rho_eff_vs_dm(conditions, cfg, font_name);

metrics_table = vertcat_nonempty({dpp_metrics, density_metrics});
prediction_table = vertcat_nonempty({dpp_predictions, density_predictions});
condition_summary = summarize_conditions(conditions);

%% Export stable manuscript figures and timestamped audit artifacts

run_timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
run_dir = fullfile(cfg.outputs.root, 'runs', run_timestamp);
if ~isfolder(run_dir)
    mkdir(run_dir);
end

if cfg.outputs.export
    if ~isfolder(cfg.outputs.root)
        mkdir(cfg.outputs.root);
    end
    export_validation_figure(figure_dpp, cfg.outputs.root, ...
        'validation_dpp_vs_da', cfg.outputs);
    export_validation_figure(figure_density, cfg.outputs.root, ...
        'validation_rho_eff_vs_dm', cfg.outputs);
end

writetable(metrics_table, fullfile(run_dir, 'validation_metrics.csv'));
writetable(prediction_table, fullfile(run_dir, 'validation_predictions.csv'));
writetable(condition_summary, fullfile(run_dir, 'condition_summary.csv'));
write_json(fullfile(run_dir, 'resolved_config.json'), cfg);

manifest = build_manifest(cfg, source_info, density_source_data.metadata, ...
    run_timestamp, font_name, condition_summary);
write_json(fullfile(run_dir, 'run_manifest.json'), manifest);

fit_curves = struct('dpp_vs_da', dpp_fits, ...
    'rho_eff_vs_dm', density_fits);
save(fullfile(run_dir, 'validation_results.mat'), 'metrics_table', ...
    'prediction_table', 'condition_summary', 'fit_curves', 'manifest', '-v7.3');

result = struct();
result.metrics = metrics_table;
result.predictions = prediction_table;
result.condition_summary = condition_summary;
result.fit_curves = fit_curves;
result.figures = struct('dpp_vs_da', figure_dpp, ...
    'rho_eff_vs_dm', figure_density);
result.run_dir = run_dir;
result.font_name = font_name;

fprintf('Validation outputs written to:\n%s\n', cfg.outputs.root);
fprintf('Run-specific metrics and provenance written to:\n%s\n', run_dir);

if strcmpi(cfg.figures.visible, 'off')
    close(figure_dpp);
    close(figure_density);
    result.figures = struct('dpp_vs_da', [], 'rho_eff_vs_dm', []);
end

end

function validate_source_payloads(ld2_data, tem_data, density_data)
%VALIDATE_SOURCE_PAYLOADS Fail early when an input has the wrong schema.

required_ld2 = {'parsdata', 'fl', 'r_n_agg'};
for i = 1:numel(required_ld2)
    if ~isfield(ld2_data, required_ld2{i})
        error('PFAL:RUN_MAIN_VALIDATION:MissingLD2Variable', ...
            'The LD2 source is missing variable "%s".', required_ld2{i});
    end
end
if ~isstruct(ld2_data.parsdata) || ...
        numel(ld2_data.parsdata) ~= numel(ld2_data.r_n_agg)
    error('PFAL:RUN_MAIN_VALIDATION:InvalidLD2Snapshots', ...
        'parsdata and r_n_agg must describe the same number of snapshots.');
end
if ~isstruct(ld2_data.fl) || ~all(isfield(ld2_data.fl, {'mu', 'lambda'}))
    error('PFAL:RUN_MAIN_VALIDATION:InvalidFluidData', ...
        'The LD2 source fluid struct must contain mu and lambda.');
end
if ~isfield(tem_data, 'aggregate_table') || ~istable(tem_data.aggregate_table)
    error('PFAL:RUN_MAIN_VALIDATION:InvalidTEMData', ...
        'The TEM source must contain aggregate_table.');
end
required_tem = {'entry_id', 'da_nm', 'dbarpp_nm'};
assert_table_variables(tem_data.aggregate_table, required_tem, 'TEM aggregate_table');
if ~isfield(density_data, 'effective_density_data') || ...
        ~istable(density_data.effective_density_data)
    error('PFAL:RUN_MAIN_VALIDATION:InvalidDensityData', ...
        'The density source must contain effective_density_data.');
end
required_density = {'condition_id', 'mobility_mode_nm', ...
    'effective_density_kg_m3'};
assert_table_variables(density_data.effective_density_data, ...
    required_density, 'effective_density_data');
if ~isfield(density_data, 'metadata') || ~isstruct(density_data.metadata)
    error('PFAL:RUN_MAIN_VALIDATION:MissingDensityMetadata', ...
        'The normalized density source must contain its metadata struct.');
end

end

function conditions = prepare_conditions(configured_conditions, ld2_data, ...
        tem_table, density_table, material_density)
%PREPARE_CONDITIONS Resolve each enabled name to simulation and experiment.

enabled = configured_conditions([configured_conditions.enabled]);
conditions = repmat(struct(), numel(enabled), 1);
for i = 1:numel(enabled)
    condition_cfg = enabled(i);
    simulation = select_simulation_population(ld2_data.parsdata, ...
        ld2_data.r_n_agg, condition_cfg.ld2_fractions);
    simulation = derive_simulation_observables(simulation, ld2_data.fl, ...
        material_density);

    conditions(i).id = condition_cfg.id;
    conditions(i).label = condition_cfg.label;
    conditions(i).include = condition_cfg.include;
    conditions(i).style = condition_cfg.style;
    conditions(i).simulation = simulation;
    conditions(i).tem = table();
    conditions(i).effective_density = table();

    if condition_cfg.include.dpp_vs_da
        mask = strcmp(string(tem_table.entry_id), ...
            string(condition_cfg.tem_entry_id));
        conditions(i).tem = tem_table(mask, :);
        if isempty(conditions(i).tem)
            error('PFAL:RUN_MAIN_VALIDATION:UnknownTEMCondition', ...
                'No TEM rows matched entry_id "%s" for condition "%s".', ...
                condition_cfg.tem_entry_id, condition_cfg.id);
        end
    end

    if condition_cfg.include.rho_eff_vs_dm
        mask = strcmp(string(density_table.condition_id), ...
            string(condition_cfg.effective_density_condition_id));
        conditions(i).effective_density = density_table(mask, :);
        if isempty(conditions(i).effective_density)
            error('PFAL:RUN_MAIN_VALIDATION:UnknownDensityCondition', ...
                ['No effective-density rows matched condition_id "%s" ' ...
                 'for condition "%s".'], ...
                condition_cfg.effective_density_condition_id, condition_cfg.id);
        end
    end
end

end

function population = select_simulation_population(parsdata, fractions, requested)
%SELECT_SIMULATION_POPULATION Match fractions without relying on row numbers.

indices = zeros(numel(requested), 1);
fractions = double(fractions(:));
for i = 1:numel(requested)
    target = requested(i);
    tolerance = max(1e-12, 1e-9 * abs(target));
    matched = find(abs(fractions - target) <= tolerance);
    if numel(matched) ~= 1
        error('PFAL:RUN_MAIN_VALIDATION:InvalidLD2FractionMapping', ...
            ['Expected exactly one LD2 snapshot for fraction %.12g but ' ...
             'found %d. Available fractions: %s'], ...
            target, numel(matched), mat2str(fractions.'));
    end
    indices(i) = matched;
end

selected = parsdata(indices);
while numel(selected) > 1
    % The established merge helper compares complete primary-ID sets before
    % removing later duplicates from nested LD2 snapshots.
    selected = UTILS.MERGE_PARSDATA_ROWS(selected, 1);
end
population = selected(1);

required_fields = {'dpp', 'da', 'dg', 'pp'};
for i = 1:numel(required_fields)
    if ~isfield(population, required_fields{i})
        error('PFAL:RUN_MAIN_VALIDATION:MissingPopulationField', ...
            'The selected LD2 population is missing field "%s".', ...
            required_fields{i});
    end
end

end

function simulation = derive_simulation_observables(population, fl, material_density)
%DERIVE_SIMULATION_OBSERVABLES Recalculate mobility diameter and density.

simulation = struct();
simulation.da_nm = 1e9 * double(population.da(:));
simulation.dpp_nm = 1e9 * double(population.dpp(:));

% Later LD2 snapshots saved da into dm, so dm is deliberately never read.
dm_m = TRANSP.DIAMOBIL(double(population.dg(:)), ...
    double(population.da(:)), fl);
simulation.dm_nm = 1e9 * dm_m(:);

particle_volume_factor = cellfun(@(particles) ...
    (pi / 6) * sum(double(particles(:, 2)).^3), population.pp(:));
simulation.rho_eff_kg_m3 = material_density .* ...
    particle_volume_factor(:) ./ (dm_m(:).^3);
simulation.aggregate_count = numel(simulation.da_nm);

assert_positive_finite(simulation.da_nm, 'simulated aerodynamic diameter');
assert_positive_finite(simulation.dpp_nm, 'simulated primary-particle diameter');
assert_positive_finite(simulation.dm_nm, 'recalculated mobility diameter');
assert_positive_finite(simulation.rho_eff_kg_m3, 'simulated effective density');

end

function [fig, metrics, predictions, fit_curves] = plot_dpp_vs_da( ...
        conditions, cfg, font_name)
%PLOT_DPP_VS_DA Build the primary-particle morphology validation figure.

panel_conditions = conditions(arrayfun(@(x) x.include.dpp_vs_da, conditions));
if isempty(panel_conditions)
    error('PFAL:RUN_MAIN_VALIDATION:NoDppConditions', ...
        'No enabled conditions participate in dpp_vs_da.');
end

fig = create_figure(cfg.figures, cfg.figures.dpp_vs_da, font_name);
ax = axes(fig);
hold(ax, 'on');
configure_axes(ax, font_name, cfg.figures.font);
axis_text = scientific_axis_text(cfg.figures.font.interpreter);
xlabel(ax, axis_text.da, 'Interpreter', cfg.figures.font.interpreter, ...
    'FontName', font_name, 'FontSize', cfg.figures.font.label_size);
ylabel(ax, axis_text.dpp, 'Interpreter', cfg.figures.font.interpreter, ...
    'FontName', font_name, 'FontSize', cfg.figures.font.label_size);

fit_handles = gobjects(0);
experiment_handles = gobjects(0);
fit_labels = {};
experiment_labels = {};
metric_cells = cell(numel(panel_conditions), 1);
prediction_cells = cell(numel(panel_conditions), 1);
fit_curves = repmat(empty_fit_record(), numel(panel_conditions), 1);

rng(cfg.fits.dpp_vs_da.rng_seed, 'twister');
for i = 1:numel(panel_conditions)
    condition = panel_conditions(i);
    color = hex_color(condition.style.color);
    simulation = condition.simulation;

    if cfg.figures.show_raw_simulation
        scatter(ax, simulation.da_nm, simulation.dpp_nm, ...
            condition.style.simulation_marker_size, color, ...
            condition.style.simulation_marker, 'filled', ...
            'MarkerFaceAlpha', cfg.figures.simulation_marker_alpha, ...
            'MarkerEdgeAlpha', cfg.figures.simulation_marker_alpha, ...
            'HandleVisibility', 'off');
    end

    fit_record = calculate_fit(simulation.da_nm, simulation.dpp_nm, ...
        cfg.fits.dpp_vs_da, condition.id);
    fit_curves(i) = fit_record;
    [fit_handles(end + 1), fit_labels{end + 1}] = plot_fit_or_raw( ...
        ax, fit_record, simulation.da_nm, simulation.dpp_nm, ...
        condition, color, cfg.figures); %#ok<AGROW>

    experiment = condition.tem;
    if cfg.figures.show_tem_error_bars && ...
            all(ismember({'dbarpp_ci95_low_nm', 'dbarpp_ci95_high_nm'}, ...
            experiment.Properties.VariableNames))
        lower_error = experiment.dbarpp_nm - experiment.dbarpp_ci95_low_nm;
        upper_error = experiment.dbarpp_ci95_high_nm - experiment.dbarpp_nm;
        experiment_handles(end + 1) = errorbar(ax, experiment.da_nm, ...
            experiment.dbarpp_nm, lower_error, upper_error, ...
            'LineStyle', 'none', 'Color', color, ...
            'Marker', condition.style.experimental_marker, ...
            'MarkerSize', sqrt(condition.style.experimental_marker_size), ...
            'LineWidth', 1.1); %#ok<AGROW>
    else
        experiment_handles(end + 1) = scatter(ax, experiment.da_nm, ...
            experiment.dbarpp_nm, condition.style.experimental_marker_size, ...
            color, condition.style.experimental_marker, ...
            'LineWidth', 1.2); %#ok<AGROW>
    end
    experiment_labels{end + 1} = sprintf('TEM: %s', condition.label); %#ok<AGROW>

    [metric_cells{i}, prediction_cells{i}] = evaluate_experiment( ...
        condition.id, 'dpp_vs_da', experiment.da_nm, ...
        experiment.dbarpp_nm, fit_record);
end

[reference_handle, reference_label] = plot_dpp_reference(ax, ...
    panel_conditions, cfg.references.dpp_vs_da);
apply_limits(ax, panel_conditions, 'dpp_vs_da', cfg.figures.dpp_vs_da);
make_legend(ax, fit_handles, experiment_handles, reference_handle, ...
    fit_labels, experiment_labels, reference_label, cfg.figures.font, font_name);

metrics = vertcat(metric_cells{:});
predictions = vertcat(prediction_cells{:});

end

function [fig, metrics, predictions, fit_curves] = plot_rho_eff_vs_dm( ...
        conditions, cfg, font_name)
%PLOT_RHO_EFF_VS_DM Build the effective-density validation figure.

panel_conditions = conditions(arrayfun(@(x) x.include.rho_eff_vs_dm, conditions));
if isempty(panel_conditions)
    error('PFAL:RUN_MAIN_VALIDATION:NoDensityConditions', ...
        'No enabled conditions participate in rho_eff_vs_dm.');
end

fig = create_figure(cfg.figures, cfg.figures.rho_eff_vs_dm, font_name);
ax = axes(fig);
hold(ax, 'on');
configure_axes(ax, font_name, cfg.figures.font);
axis_text = scientific_axis_text(cfg.figures.font.interpreter);
xlabel(ax, axis_text.dm, 'Interpreter', cfg.figures.font.interpreter, ...
    'FontName', font_name, 'FontSize', cfg.figures.font.label_size);
ylabel(ax, axis_text.rho_eff, ...
    'Interpreter', cfg.figures.font.interpreter, ...
    'FontName', font_name, 'FontSize', cfg.figures.font.label_size);

fit_handles = gobjects(0);
experiment_handles = gobjects(0);
fit_labels = {};
experiment_labels = {};
metric_cells = cell(numel(panel_conditions), 1);
prediction_cells = cell(numel(panel_conditions), 1);
fit_curves = repmat(empty_fit_record(), numel(panel_conditions), 1);

rng(cfg.fits.rho_eff_vs_dm.rng_seed, 'twister');
for i = 1:numel(panel_conditions)
    condition = panel_conditions(i);
    color = hex_color(condition.style.color);
    simulation = condition.simulation;

    if cfg.figures.show_raw_simulation
        scatter(ax, simulation.dm_nm, simulation.rho_eff_kg_m3, ...
            condition.style.simulation_marker_size, color, ...
            condition.style.simulation_marker, 'filled', ...
            'MarkerFaceAlpha', cfg.figures.simulation_marker_alpha, ...
            'MarkerEdgeAlpha', cfg.figures.simulation_marker_alpha, ...
            'HandleVisibility', 'off');
    end

    fit_record = calculate_fit(simulation.dm_nm, ...
        simulation.rho_eff_kg_m3, cfg.fits.rho_eff_vs_dm, condition.id);
    fit_curves(i) = fit_record;
    [fit_handles(end + 1), fit_labels{end + 1}] = plot_fit_or_raw( ...
        ax, fit_record, simulation.dm_nm, simulation.rho_eff_kg_m3, ...
        condition, color, cfg.figures); %#ok<AGROW>

    experiment = condition.effective_density;
    experiment_handles(end + 1) = scatter(ax, ...
        experiment.mobility_mode_nm, experiment.effective_density_kg_m3, ...
        condition.style.experimental_marker_size, color, ...
        condition.style.experimental_marker, 'LineWidth', 1.2); %#ok<AGROW>
    experiment_labels{end + 1} = sprintf('AAC-SMPS: %s', condition.label); %#ok<AGROW>

    [metric_cells{i}, prediction_cells{i}] = evaluate_experiment( ...
        condition.id, 'rho_eff_vs_dm', experiment.mobility_mode_nm, ...
        experiment.effective_density_kg_m3, fit_record);
end

[reference_handle, reference_label] = plot_density_reference(ax, ...
    panel_conditions, cfg.references.rho_eff_vs_dm);
apply_limits(ax, panel_conditions, 'rho_eff_vs_dm', ...
    cfg.figures.rho_eff_vs_dm);
make_legend(ax, fit_handles, experiment_handles, reference_handle, ...
    fit_labels, experiment_labels, reference_label, cfg.figures.font, font_name);

metrics = vertcat(metric_cells{:});
predictions = vertcat(prediction_cells{:});

end

function fit_record = calculate_fit(x, y, fit_cfg, condition_id)
%CALCULATE_FIT Run the configured Bayesian model or record its disabled state.

fit_record = empty_fit_record();
fit_record.condition_id = condition_id;
fit_record.enabled = fit_cfg.enabled;
fit_record.show_band = fit_cfg.show_band;
fit_record.x_transform = fit_cfg.x_transform;
fit_record.credible_level = fit_cfg.credible_level;
fit_record.x_support = [min(x), max(x)];
if ~fit_cfg.enabled
    fit_record.status = 'fit_disabled';
    return
end

[fit_record.y, fit_record.x, fit_record.bounds_y, ...
    fit_record.slope, fit_record.bounds_slope, diagnostics] = ...
    UTILS.FIT_POLY(x, y, ones(numel(x), 1), fit_cfg.resolution, ...
    'Degree', fit_cfg.degree, ...
    'XTransform', fit_cfg.x_transform, ...
    'YTransform', fit_cfg.y_transform, ...
    'NSamples', fit_cfg.posterior_samples, ...
    'SlopeOffset', fit_cfg.slope_offset, ...
    'CredibleLevel', fit_cfg.credible_level, ...
    'PriorMu', fit_cfg.prior_mu, ...
    'PriorVScale', fit_cfg.prior_v_scale, ...
    'PriorA', fit_cfg.prior_a, ...
    'PriorB', fit_cfg.prior_b);
fit_record.status = 'fit_enabled';
fit_record.posterior = diagnostics.Posterior;
fit_record.options = diagnostics.Options;

end

function record = empty_fit_record()
%EMPTY_FIT_RECORD Provide a consistent serializable fit result shape.

record = struct('condition_id', '', 'enabled', false, 'show_band', false, ...
    'status', 'fit_disabled', 'x_transform', 'log10', ...
    'credible_level', NaN, 'x_support', [NaN, NaN], ...
    'x', [], 'y', [], 'bounds_y', [], 'slope', [], ...
    'bounds_slope', [], 'posterior', struct(), 'options', struct());

end

function [handle, label] = plot_fit_or_raw(ax, fit_record, x, y, ...
        condition, color, figure_cfg)
%PLOT_FIT_OR_RAW Use the raw cloud as the legend item when fitting is off.

if fit_record.enabled
    if fit_record.show_band && figure_cfg.band_alpha > 0 && ...
            ~isempty(fit_record.bounds_y)
        fill(ax, [fit_record.x; flipud(fit_record.x)], ...
            [fit_record.bounds_y(:, 1); flipud(fit_record.bounds_y(:, 2))], ...
            color, 'EdgeColor', 'none', ...
            'FaceAlpha', figure_cfg.band_alpha, 'HandleVisibility', 'off');
    end
    handle = plot(ax, fit_record.x, fit_record.y, 'Color', color, ...
        'LineWidth', condition.style.fit_line_width);
    label = sprintf('LD2 fit: %s', condition.label);
else
    handle = scatter(ax, x, y, condition.style.simulation_marker_size, ...
        color, condition.style.simulation_marker, 'filled', ...
        'MarkerFaceAlpha', max(figure_cfg.simulation_marker_alpha, 0.35), ...
        'MarkerEdgeAlpha', max(figure_cfg.simulation_marker_alpha, 0.35));
    label = sprintf('LD2: %s (fit off)', condition.label);
end

end

function [metric_row, prediction_rows] = evaluate_experiment( ...
        condition_id, panel_id, experimental_x, experimental_y, fit_record)
%EVALUATE_EXPERIMENT Compare observations only within simulation support.

experimental_x = double(experimental_x(:));
experimental_y = double(experimental_y(:));
n_points = numel(experimental_x);
predicted = nan(n_points, 1);
lower = nan(n_points, 1);
upper = nan(n_points, 1);
residual = nan(n_points, 1);
within_band = nan(n_points, 1);
in_support = experimental_x >= fit_record.x_support(1) & ...
    experimental_x <= fit_record.x_support(2);

if fit_record.enabled
    predicted(in_support) = interpolate_curve(fit_record.x, fit_record.y, ...
        experimental_x(in_support), fit_record.x_transform);
    lower(in_support) = interpolate_curve(fit_record.x, ...
        fit_record.bounds_y(:, 1), experimental_x(in_support), ...
        fit_record.x_transform);
    upper(in_support) = interpolate_curve(fit_record.x, ...
        fit_record.bounds_y(:, 2), experimental_x(in_support), ...
        fit_record.x_transform);
    valid = in_support & experimental_y > 0 & predicted > 0;
    residual(valid) = log10(experimental_y(valid)) - log10(predicted(valid));
    within_band(valid) = double(experimental_y(valid) >= lower(valid) & ...
        experimental_y(valid) <= upper(valid));
else
    valid = false(n_points, 1);
end

if any(valid)
    log_bias = mean(residual(valid));
    log_mae = mean(abs(residual(valid)));
    log_rmse = sqrt(mean(residual(valid).^2));
    band_coverage = mean(within_band(valid));
else
    log_bias = NaN;
    log_mae = NaN;
    log_rmse = NaN;
    band_coverage = NaN;
end

metric_row = table(string(condition_id), string(panel_id), ...
    string(fit_record.status), n_points, sum(valid), sum(~in_support), ...
    log_bias, log_mae, log_rmse, band_coverage, ...
    'VariableNames', {'condition_id', 'panel_id', 'fit_status', ...
    'n_experimental', 'n_in_support', 'n_outside_support', ...
    'log10_bias', 'log10_mae', 'log10_rmse', ...
    'credible_band_coverage'});

prediction_rows = table( ...
    repmat(string(condition_id), n_points, 1), ...
    repmat(string(panel_id), n_points, 1), ...
    (1:n_points)', experimental_x, experimental_y, predicted, lower, upper, ...
    in_support, residual, within_band, ...
    'VariableNames', {'condition_id', 'panel_id', ...
    'experimental_point_index', 'x_value', 'observed_y', 'predicted_y', ...
    'credible_lower', 'credible_upper', 'in_simulation_support', ...
    'log10_residual', 'within_credible_band'});

end

function values = interpolate_curve(curve_x, curve_y, query_x, x_transform)
%INTERPOLATE_CURVE Evaluate predictions on the configured predictor scale.

if strcmpi(x_transform, 'log10')
    values = interp1(log10(curve_x), curve_y, log10(query_x), 'linear', NaN);
else
    values = interp1(curve_x, curve_y, query_x, 'linear', NaN);
end

end

function [handle, label] = plot_dpp_reference(ax, conditions, reference)
%PLOT_DPP_REFERENCE Draw the configurable morphology literature relation.

handle = gobjects(0);
label = {};
if ~reference.enabled
    return
end
x = collect_panel_values(conditions, 'dpp_vs_da', 'x');
x_reference = logspace(log10(min(x) / 1.05), log10(max(x) * 1.05), 200)';
y_reference = reference.value_at_100_nm .* ...
    (x_reference ./ 100).^reference.exponent;
handle = plot(ax, x_reference, y_reference, ...
    'Color', hex_color(reference.color), ...
    'LineStyle', reference.line_style, 'LineWidth', reference.line_width);
label = {reference.label};

end

function [handle, label] = plot_density_reference(ax, conditions, reference)
%PLOT_DENSITY_REFERENCE Draw the configurable mass-mobility relation.

handle = gobjects(0);
label = {};
if ~reference.enabled
    return
end
x = collect_panel_values(conditions, 'rho_eff_vs_dm', 'x');
x_reference = logspace(log10(min(x) / 1.05), log10(max(x) * 1.05), 200)';
y_reference = reference.value_at_100_nm .* ...
    (x_reference ./ 100).^(reference.mass_mobility_exponent - 3);
handle = plot(ax, x_reference, y_reference, ...
    'Color', hex_color(reference.color), ...
    'LineStyle', reference.line_style, 'LineWidth', reference.line_width);
label = {reference.label};

end

function values = collect_panel_values(conditions, panel_id, axis_id)
%COLLECT_PANEL_VALUES Combine simulation and experiment for automatic limits.

values = [];
for i = 1:numel(conditions)
    condition = conditions(i);
    switch panel_id
        case 'dpp_vs_da'
            if strcmp(axis_id, 'x')
                values = [values; condition.simulation.da_nm; ...
                    condition.tem.da_nm]; %#ok<AGROW>
            else
                values = [values; condition.simulation.dpp_nm; ...
                    condition.tem.dbarpp_nm]; %#ok<AGROW>
            end
        case 'rho_eff_vs_dm'
            if strcmp(axis_id, 'x')
                values = [values; condition.simulation.dm_nm; ...
                    condition.effective_density.mobility_mode_nm]; %#ok<AGROW>
            else
                values = [values; condition.simulation.rho_eff_kg_m3; ...
                    condition.effective_density.effective_density_kg_m3]; %#ok<AGROW>
            end
    end
end
values = double(values(:));
values = values(isfinite(values) & values > 0);

end

function apply_limits(ax, conditions, panel_id, panel_cfg)
%APPLY_LIMITS Use configured limits or padded positive data extents.

x = collect_panel_values(conditions, panel_id, 'x');
y = collect_panel_values(conditions, panel_id, 'y');
if isempty(panel_cfg.x_limits)
    xlim(ax, [min(x) / 1.08, max(x) * 1.08]);
else
    xlim(ax, panel_cfg.x_limits);
end
if isempty(panel_cfg.y_limits)
    ylim(ax, [min(y) / 1.08, max(y) * 1.08]);
else
    ylim(ax, panel_cfg.y_limits);
end

end

function fig = create_figure(figures_cfg, panel_cfg, font_name)
%CREATE_FIGURE Initialize a publication figure with configured typography.

fig = figure('Visible', figures_cfg.visible, 'Color', 'white', ...
    'Position', panel_cfg.position, 'DefaultTextFontName', font_name, ...
    'DefaultAxesFontName', font_name);

end

function configure_axes(ax, font_name, font_cfg)
%CONFIGURE_AXES Apply shared log-axis styling to a validation panel.

set(ax, 'XScale', 'log', 'YScale', 'log', 'Box', 'on', ...
    'TickLength', [0.02, 0.02], 'FontName', font_name, ...
    'FontSize', font_cfg.axis_size, 'FontWeight', font_cfg.weight, ...
    'TickLabelInterpreter', font_cfg.interpreter, 'Layer', 'top');

end

function labels = scientific_axis_text(interpreter)
%SCIENTIFIC_AXIS_TEXT Italicize variables while keeping subscripts upright.

switch lower(interpreter)
    case 'latex'
        labels.da = '$d_{\mathrm{a}}$ [nm]';
        labels.dpp = '$d_{\mathrm{pp}}$ [nm]';
        labels.dm = '$d_{\mathrm{m}}$ [nm]';
        labels.rho_eff = '$\rho_{\mathrm{eff}}$ [kg m$^{-3}$]';
    case 'tex'
        labels.da = '{\it d}_{\rm a} [nm]';
        labels.dpp = '{\it d}_{\rm pp} [nm]';
        labels.dm = '{\it d}_{\rm m} [nm]';
        labels.rho_eff = '{\it \rho}_{\rm eff} [kg m^{-3}]';
    otherwise
        labels.da = 'd_a [nm]';
        labels.dpp = 'd_pp [nm]';
        labels.dm = 'd_m [nm]';
        labels.rho_eff = 'rho_eff [kg m^-3]';
end

end

function make_legend(ax, fit_handles, experiment_handles, reference_handle, ...
        fit_labels, experiment_labels, reference_label, font_cfg, font_name)
%MAKE_LEGEND Assemble simulation, experiment, and reference entries.

handles = [fit_handles(:); experiment_handles(:); reference_handle(:)];
labels = [fit_labels(:); experiment_labels(:); reference_label(:)];
legend(ax, handles, labels, 'Interpreter', font_cfg.interpreter, ...
    'FontName', font_name, 'FontSize', font_cfg.legend_size, ...
    'FontWeight', font_cfg.weight, 'Location', 'northoutside', ...
    'NumColumns', 2, 'Orientation', 'horizontal');

end

function font_name = resolve_font(font_cfg)
%RESOLVE_FONT Use Segoe UI when installed and warn before falling back.

available_fonts = listfonts;
matched = find(strcmpi(available_fonts, font_cfg.family), 1);
if ~isempty(matched)
    font_name = available_fonts{matched};
    return
end

fallback = find(strcmpi(available_fonts, font_cfg.fallback), 1);
if ~isempty(fallback)
    font_name = available_fonts{fallback};
else
    font_name = get(groot, 'DefaultAxesFontName');
end
warning('PFAL:RUN_MAIN_VALIDATION:FontFallback', ...
    'Font "%s" is unavailable; using "%s".', font_cfg.family, font_name);

end

function color = hex_color(value)
%HEX_COLOR Convert a six-digit hexadecimal color to MATLAB RGB values.

value = char(value);
if numel(value) == 7 && value(1) == '#'
    value = value(2:end);
end
if numel(value) ~= 6 || any(~ismember(lower(value), '0123456789abcdef'))
    error('PFAL:RUN_MAIN_VALIDATION:InvalidColor', ...
        'Expected a six-digit hexadecimal color, received "%s".', value);
end
color = reshape(sscanf(value, '%2x'), 1, 3) ./ 255;

end

function summary = summarize_conditions(conditions)
%SUMMARIZE_CONDITIONS Record exact simulation and experimental sample counts.

n = numel(conditions);
condition_id = strings(n, 1);
condition_label = strings(n, 1);
n_simulation = zeros(n, 1);
n_tem = zeros(n, 1);
n_effective_density = zeros(n, 1);
for i = 1:n
    condition_id(i) = string(conditions(i).id);
    condition_label(i) = string(conditions(i).label);
    n_simulation(i) = conditions(i).simulation.aggregate_count;
    n_tem(i) = height(conditions(i).tem);
    n_effective_density(i) = height(conditions(i).effective_density);
end
summary = table(condition_id, condition_label, n_simulation, n_tem, ...
    n_effective_density);

end

function manifest = build_manifest(cfg, source_info, density_metadata, ...
        run_timestamp, font_name, condition_summary)
%BUILD_MANIFEST Capture source identity and resolved run behavior.

manifest = struct();
manifest.schema_version = '1.0.0';
manifest.generated_at = char(datetime('now', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ss'));
manifest.run_timestamp = run_timestamp;
manifest.config_file = cfg.config_file;
manifest.font_requested = cfg.figures.font.family;
manifest.font_used = font_name;
manifest.mobility_diameter_source = ...
    'Recalculated with TRANSP.DIAMOBIL from parsdata.dg, parsdata.da, and fl.';
manifest.sources = struct( ...
    'ld2', source_manifest(source_info.ld2), ...
    'tem', source_manifest(source_info.tem), ...
    'effective_density', source_manifest(source_info.effective_density));
manifest.effective_density_source_metadata = density_metadata;
manifest.condition_counts = table2struct(condition_summary);

end

function output = source_manifest(source)
%SOURCE_MANIFEST Add file size, timestamp, and digest to source config data.

file_info = dir(source.resolved_file);
output = struct('id', source.id, 'resolved_file', source.resolved_file, ...
    'url', source.url, 'doi', source.doi, ...
    'sha256', UTILS.FILE_SHA256(source.resolved_file), ...
    'bytes', file_info.bytes, ...
    'last_modified', char(datetime(file_info.datenum, ...
    'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd''T''HH:mm:ss')), ...
    'loaded_variables', {source.loaded_variables});

end

function export_validation_figure(fig, output_root, base_name, output_cfg)
%EXPORT_VALIDATION_FIGURE Write stable PDF, PNG, and optional FIG artifacts.

if output_cfg.pdf
    exportgraphics(fig, fullfile(output_root, [base_name, '.pdf']), ...
        'ContentType', 'vector');
end
if output_cfg.png
    exportgraphics(fig, fullfile(output_root, [base_name, '.png']), ...
        'Resolution', output_cfg.png_resolution);
end
if output_cfg.save_figures
    savefig(fig, fullfile(output_root, [base_name, '.fig']));
end

end

function write_json(file_path, value)
%WRITE_JSON Write readable JSON while supporting older MATLAB releases.

try
    text = jsonencode(value, 'PrettyPrint', true);
catch
    text = jsonencode(value);
end
fid = fopen(file_path, 'w');
if fid < 0
    error('PFAL:RUN_MAIN_VALIDATION:JSONWriteFailed', ...
        'Could not open JSON output for writing: %s', file_path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', text);

end

function assert_table_variables(input_table, required, table_label)
%ASSERT_TABLE_VARIABLES Validate table schema with an actionable error.

missing = setdiff(required, input_table.Properties.VariableNames);
if ~isempty(missing)
    error('PFAL:RUN_MAIN_VALIDATION:MissingTableVariables', ...
        '%s is missing required variables: %s', ...
        table_label, strjoin(missing, ', '));
end

end

function assert_positive_finite(values, label)
%ASSERT_POSITIVE_FINITE Protect log plots and log-space validation metrics.

if isempty(values) || any(~isfinite(values)) || any(values <= 0)
    error('PFAL:RUN_MAIN_VALIDATION:InvalidObservable', ...
        'The %s values must be non-empty, finite, and positive.', label);
end

end

function table_out = vertcat_nonempty(table_cells)
%VERTCAT_NONEMPTY Combine table results while tolerating an empty panel result.

mask = cellfun(@(value) istable(value) && ~isempty(value), table_cells);
if any(mask)
    table_out = vertcat(table_cells{mask});
else
    table_out = table();
end

end
