function [effective_density_data, metadata] = NORMALIZE_EFFECTIVE_DENSITY_DATA( ...
        source_file, output_file, varargin)
%NORMALIZE_EFFECTIVE_DENSITY_DATA Create a compact, lossless validation MAT.
%   This utility copies the grouped ODIAS results into a long-format table.
%   It does not rerun the tandem inversion, mode selection, density equation,
%   or any fit. The original MAT remains the numerical source of record.
%
%   Optional name-value input:
%     'ExpectedSourceSHA256' - reject an unexpected source artifact.

expected_hash = '';
if mod(numel(varargin), 2) ~= 0
    error('PFAL:NORMALIZE_EFFECTIVE_DENSITY_DATA:BadOptions', ...
        'Optional inputs must be supplied as name-value pairs.');
end
for i = 1:2:numel(varargin)
    switch lower(char(varargin{i}))
        case 'expectedsourcesha256'
            expected_hash = lower(char(varargin{i + 1}));
        otherwise
            error('PFAL:NORMALIZE_EFFECTIVE_DENSITY_DATA:UnknownOption', ...
                'Unknown option: %s', char(varargin{i}));
    end
end

source_file = char(source_file);
output_file = char(output_file);
if ~isfile(source_file)
    error('PFAL:NORMALIZE_EFFECTIVE_DENSITY_DATA:MissingSource', ...
        'Effective-density source MAT was not found: %s', source_file);
end

source_hash_before = UTILS.FILE_SHA256(source_file);
if ~isempty(expected_hash) && ~strcmpi(source_hash_before, expected_hash)
    error('PFAL:NORMALIZE_EFFECTIVE_DENSITY_DATA:SourceChecksumMismatch', ...
        ['The source MAT does not match the expected SHA-256.\n' ...
         'Expected: %s\nActual:   %s'], expected_hash, source_hash_before);
end

% Loading named variables avoids bringing the large raw distribution arrays
% and saved graphics objects into the validation preparation step.
source = load(source_file, 'dist_grp', 'test_condition');
if ~isfield(source, 'dist_grp') || ~isstruct(source.dist_grp)
    error('PFAL:NORMALIZE_EFFECTIVE_DENSITY_DATA:MissingLegacyData', ...
        'The source MAT does not contain the expected dist_grp struct.');
end

legacy_groups = source.dist_grp(:);
condition_ids = ["low_agglomeration"; "moderate_collapse"; ...
    "high_agglomeration"; "extensive_collapse"];
default_labels = ["Low agglomeration"; "Moderate collapse"; ...
    "High agglomeration"; "Extensive collapse"];
if numel(legacy_groups) ~= numel(condition_ids)
    error('PFAL:NORMALIZE_EFFECTIVE_DENSITY_DATA:UnexpectedGroupCount', ...
        'Expected four ODIAS condition groups but found %d.', numel(legacy_groups));
end

condition_labels = extract_condition_labels(source, default_labels);
required_fields = {'da', 'd_mode', 'd_gm', 'sigma_g', 'rho_eff'};
rows = cell(numel(legacy_groups), 1);

for group_index = 1:numel(legacy_groups)
    group = legacy_groups(group_index);
    for field_index = 1:numel(required_fields)
        if ~isfield(group, required_fields{field_index})
            error('PFAL:NORMALIZE_EFFECTIVE_DENSITY_DATA:MissingField', ...
                'dist_grp(%d) is missing field "%s".', ...
                group_index, required_fields{field_index});
        end
    end

    n_points = numel(group.da);
    lengths = cellfun(@(name) numel(group.(name)), required_fields);
    if any(lengths ~= n_points)
        error('PFAL:NORMALIZE_EFFECTIVE_DENSITY_DATA:LengthMismatch', ...
            'The numeric fields in dist_grp(%d) do not have equal lengths.', ...
            group_index);
    end

    rows{group_index} = table( ...
        repmat(condition_ids(group_index), n_points, 1), ...
        repmat(condition_labels(group_index), n_points, 1), ...
        repmat(double(group_index), n_points, 1), ...
        (1:n_points)', ...
        double(group.da(:)), ...
        double(group.d_mode(:)), ...
        double(group.d_gm(:)), ...
        double(group.sigma_g(:)), ...
        double(group.rho_eff(:)), ...
        'VariableNames', {'condition_id', 'condition_label', ...
        'source_group_index', 'point_index_within_condition', ...
        'aerodynamic_setpoint_nm', 'mobility_mode_nm', ...
        'mobility_geometric_mean_nm', ...
        'mobility_geometric_standard_deviation', ...
        'effective_density_kg_m3'});
end

effective_density_data = vertcat(rows{:});
verify_round_trip(effective_density_data, legacy_groups);

metadata = struct();
metadata.schema_version = '1.0.0';
metadata.generated_at = char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss'));
metadata.source_file = source_file;
metadata.source_sha256 = source_hash_before;
metadata.source_variable = 'dist_grp';
metadata.source_repository = 'https://github.com/Hamed-NKR/odias';
metadata.source_branch = 'developer-HN';
metadata.source_commit = 'c654348abf11599e24b63374622c31d8458cd7e4';
metadata.source_script = 'main_dma_HN_2d_v2.m';
metadata.processing_note = [ ...
    'Values were copied from dist_grp without recalculation, rounding, ' ...
    'unit conversion, inversion, mode selection, or fitting.'];
metadata.group_counts = arrayfun(@(x) numel(x.da), legacy_groups).';
metadata.legacy_field_mapping = struct( ...
    'da', 'aerodynamic_setpoint_nm', ...
    'd_mode', 'mobility_mode_nm', ...
    'd_gm', 'mobility_geometric_mean_nm', ...
    'sigma_g', 'mobility_geometric_standard_deviation', ...
    'rho_eff', 'effective_density_kg_m3');
metadata.units = struct( ...
    'aerodynamic_setpoint_nm', 'nm', ...
    'mobility_mode_nm', 'nm', ...
    'mobility_geometric_mean_nm', 'nm', ...
    'mobility_geometric_standard_deviation', '1', ...
    'effective_density_kg_m3', 'kg/m^3');
metadata.numerical_verification = struct( ...
    'round_trip_equal', true, ...
    'source_hash_unchanged', false, ...
    'row_count', height(effective_density_data));

output_dir = fileparts(output_file);
if ~isempty(output_dir) && ~isfolder(output_dir)
    mkdir(output_dir);
end

temporary_file = [tempname(output_dir), '.mat'];
temporary_cleanup = onCleanup(@() delete_if_present(temporary_file));
save(temporary_file, 'effective_density_data', 'metadata', '-v7.3');

% Reloading the candidate artifact catches serialization or table-shape issues
% before the requested output path is replaced.
candidate = load(temporary_file, 'effective_density_data');
verify_round_trip(candidate.effective_density_data, legacy_groups);

source_hash_after = UTILS.FILE_SHA256(source_file);
if ~strcmpi(source_hash_before, source_hash_after)
    error('PFAL:NORMALIZE_EFFECTIVE_DENSITY_DATA:SourceChanged', ...
        'The source MAT changed while the companion file was being prepared.');
end
metadata.numerical_verification.source_hash_unchanged = true;
save(temporary_file, 'effective_density_data', 'metadata', '-v7.3');

[moved, message] = movefile(temporary_file, output_file, 'f');
if ~moved
    error('PFAL:NORMALIZE_EFFECTIVE_DENSITY_DATA:MoveFailed', ...
        'Could not place the normalized MAT at %s: %s', output_file, message);
end
clear temporary_cleanup

output_hash = UTILS.FILE_SHA256(output_file);
write_hash_sidecar(output_file, output_hash);
fprintf('Normalized %d effective-density observations to:\n%s\n', ...
    height(effective_density_data), output_file);
fprintf('Source SHA-256: %s\nOutput SHA-256: %s\n', ...
    source_hash_before, output_hash);

end

function labels = extract_condition_labels(source, defaults)
%EXTRACT_CONDITION_LABELS Preserve source labels when their shape is usable.

labels = defaults;
if ~isfield(source, 'test_condition') || isempty(source.test_condition)
    return
end

raw = source.test_condition;
try
    if iscell(raw)
        candidate = string(raw(:));
    elseif ischar(raw)
        candidate = string(cellstr(raw));
    elseif isstring(raw) || iscategorical(raw)
        candidate = string(raw(:));
    else
        return
    end
catch
    return
end

if numel(candidate) == numel(defaults) && all(strlength(candidate) > 0)
    labels = candidate;
end

end

function verify_round_trip(tidy_data, legacy_groups)
%VERIFY_ROUND_TRIP Require exact reconstruction of every legacy numeric array.

mapping = { ...
    'da', 'aerodynamic_setpoint_nm'; ...
    'd_mode', 'mobility_mode_nm'; ...
    'd_gm', 'mobility_geometric_mean_nm'; ...
    'sigma_g', 'mobility_geometric_standard_deviation'; ...
    'rho_eff', 'effective_density_kg_m3'};

for group_index = 1:numel(legacy_groups)
    group_rows = tidy_data(tidy_data.source_group_index == group_index, :);
    [~, order] = sort(group_rows.point_index_within_condition);
    group_rows = group_rows(order, :);
    for field_index = 1:size(mapping, 1)
        legacy_name = mapping{field_index, 1};
        tidy_name = mapping{field_index, 2};
        original = legacy_groups(group_index).(legacy_name);
        reconstructed = reshape(group_rows.(tidy_name), size(original));
        if ~strcmp(class(original), class(reconstructed)) || ...
                ~isequaln(original, reconstructed)
            error('PFAL:NORMALIZE_EFFECTIVE_DENSITY_DATA:RoundTripFailed', ...
                'Exact round-trip verification failed for dist_grp(%d).%s.', ...
                group_index, legacy_name);
        end
    end
end

end

function write_hash_sidecar(output_file, hash_text)
%WRITE_HASH_SIDECAR Store the companion checksum beside the MAT artifact.

sidecar = [output_file, '.sha256'];
fid = fopen(sidecar, 'w');
if fid < 0
    error('PFAL:NORMALIZE_EFFECTIVE_DENSITY_DATA:SidecarWriteFailed', ...
        'Could not write checksum sidecar: %s', sidecar);
end
cleanup = onCleanup(@() fclose(fid));
[~, file_name, extension] = fileparts(output_file);
fprintf(fid, '%s  %s%s\n', hash_text, file_name, extension);

end

function delete_if_present(file_path)
%DELETE_IF_PRESENT Remove an incomplete temporary artifact after an error.

if isfile(file_path)
    delete(file_path);
end

end
