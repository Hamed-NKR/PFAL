function fixture_files = CREATE_MAIN_VALID_WEBTEST_FIXTURES(output_dir)
%CREATE_MAIN_VALID_WEBTEST_FIXTURES Write compact synthetic validation inputs.
%   These fixtures exercise config-based loading, condition mapping, LD2
%   pooling, plotting, and metrics without publishing or downloading the
%   full experimental and simulation artifacts.

if nargin < 1 || isempty(output_dir)
    output_dir = fullfile(repo_root(), 'data', 'main_valid_webtest', 'server');
end
if ~isfolder(output_dir)
    mkdir(output_dir);
end

%% Synthetic LD2 snapshots

r_n_agg = [1.0, 0.4, 0.2, 0.1, 0.05];
snapshot_signatures = { ...
    {[1 2], [3 4], [5 6], [7 8], [9 10], [11 12]}; ...
    {[21 22], [23 24], [25 26], [27 28]}; ...
    {[31 32], [33 34], [35 36]}; ...
    {[101 102], [103 104], [105 106]}; ...
    {[101 102], [107 108]}};
parsdata = repmat(empty_snapshot(), numel(r_n_agg), 1);
for snapshot_index = 1:numel(r_n_agg)
    parsdata(snapshot_index) = make_snapshot( ...
        snapshot_signatures{snapshot_index}, snapshot_index);
end
fl = struct('mu', 1.84e-5, 'lambda', 66e-9);
ld2_file = fullfile(output_dir, 'ld2_fixture.mat');
save(ld2_file, 'parsdata', 'fl', 'r_n_agg');

%% Synthetic processed TEM rows

low_da = [72; 91; 126; 174];
low_dpp = 17.8 .* (low_da ./ 100).^0.35 .* [0.98; 1.02; 1.01; 0.97];
high_da = [95; 142; 205];
high_dpp = 15.8 .* (high_da ./ 100).^0.27 .* [1.02; 0.96; 1.03];
aggregate_table = table( ...
    [repmat({"2024_08_19_low_agglomeration"}, 4, 1); ...
     repmat({"2025_02_06_high_agglomeration"}, 3, 1)], ...
    [low_da; high_da], [low_dpp; high_dpp], ...
    0.94 .* [low_dpp; high_dpp], 1.06 .* [low_dpp; high_dpp], ...
    'VariableNames', {'entry_id', 'da_nm', 'dbarpp_nm', ...
    'dbarpp_ci95_low_nm', 'dbarpp_ci95_high_nm'});
tem_file = fullfile(output_dir, 'tem_fixture.mat');
save(tem_file, 'aggregate_table');

%% Synthetic normalized effective-density rows

condition_id = [repmat("low_agglomeration", 4, 1); ...
    repmat("moderate_collapse", 2, 1); ...
    repmat("high_agglomeration", 4, 1); ...
    repmat("extensive_collapse", 2, 1)];
mobility_mode_nm = [70; 100; 150; 220; 100; 160; 85; 130; 190; 270; 120; 180];
effective_density_kg_m3 = [610; 510; 410; 330; 590; 520; ...
    480; 340; 245; 185; 570; 540];
effective_density_data = table(condition_id, mobility_mode_nm, ...
    effective_density_kg_m3);
metadata = struct('schema_version', 'webtest-1.0.0', ...
    'processing_note', 'Compact synthetic data; not experimental observations.');
density_file = fullfile(output_dir, 'effective_density_fixture.mat');
save(density_file, 'effective_density_data', 'metadata');

fixture_files = struct('ld2', ld2_file, 'tem', tem_file, ...
    'effective_density', density_file);
fprintf('Main validation web-test fixtures written to:\n%s\n', output_dir);

end

function snapshot = empty_snapshot()
%EMPTY_SNAPSHOT Match the fields required by the established merge helper.

snapshot = struct('dpp', [], 'sigmapp', [], 'da', [], 'dm', [], ...
    'dg', [], 'n_hyb', [], 'pp', {{}}, 'npp', []);

end

function snapshot = make_snapshot(signatures, snapshot_index)
%MAKE_SNAPSHOT Construct plausible positive aggregate observables.

n = numel(signatures);
snapshot = empty_snapshot();
snapshot.pp = cell(n, 1);
snapshot.dpp = zeros(n, 1);
snapshot.sigmapp = 1.12 * ones(n, 1);
snapshot.da = zeros(n, 1);
snapshot.dm = zeros(n, 1);
snapshot.dg = zeros(n, 1);
snapshot.n_hyb = max(snapshot_index - 1, 0) * ones(n, 1);
snapshot.npp = zeros(n, 1);

for aggregate_index = 1:n
    ids = signatures{aggregate_index}(:);
    dpp_m = (15 + 0.8 * aggregate_index + 0.3 * snapshot_index) * 1e-9;
    diameters = dpp_m .* [0.96; 1.04];
    snapshot.pp{aggregate_index} = [double(ids), diameters];
    snapshot.dpp(aggregate_index) = exp(mean(log(diameters)));
    snapshot.npp(aggregate_index) = numel(ids);
    snapshot.da(aggregate_index) = ...
        (62 + 24 * aggregate_index + 8 * snapshot_index) * 1e-9;
    snapshot.dg(aggregate_index) = 1.38 * snapshot.da(aggregate_index);
    snapshot.dm(aggregate_index) = snapshot.da(aggregate_index); % Known LD2 writer behavior.
end

end

function root_dir = repo_root()
utils_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(utils_dir);
end
