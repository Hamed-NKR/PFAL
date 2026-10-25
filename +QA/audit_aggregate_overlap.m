function report = audit_aggregate_overlap(repoRoot, varargin)
% AUDIT_AGGREGATE_OVERLAP  Geometry and index-integrity audit of the
%   PFAL two-stage (scatter -> LD2) workflow.
%
%   report = QA.AUDIT_AGGREGATE_OVERLAP(repoRoot, Name, Value, ...)
%
%   Measures primary-particle overlap (internal to each aggregate cell and
%   between independent aggregate cells, with and without periodic
%   minimum-image displacement) at every saved stage:
%     1. candidate raw pre-scatter library (pp0),
%     2. post-scatter library consumed by LD2 (pars_out.pp),
%     3. parsdata(1:5) saved LD2 populations,
%     4. every LD2 checkpoint (pars_LD2.pp),
%     5. final pars_LD2.pp.
%   It also audits column-1 primary IDs, column-6 subaggregate labels,
%   aggregate-level array alignment, matches post-scatter aggregates back to
%   the candidate pre-scatter library, verifies the post-scatter ->
%   parsdata(1) transformation is rigid, localizes the earliest overlap in
%   time, and runs isolated property tests against COL.CONNECT, COL.GROW,
%   and COL.UNITE using in-memory aggregate fixtures.
%
%   Overlap detection is independent of COL.OVR. The implementation uses an
%   exact voxel-hash broad phase validated against brute force and requires
%   no additional toolboxes.
%
%   Options (defaults point at the 2026-05-16 from-main-scatter LD2 run):
%     'PreScatterFile'       data/main_scatter/LD1__gamma_1_35__merge.mat
%     'PostScatterFile'      data/main_scatter/scaled_aggs_for_LD2_from_main_scatter.mat
%     'LD2FinalFile'         results/.../LD2__from_main_scatter__final__2026-06-28_20-05-29.mat
%     'CheckpointFiles'      'auto' (all *checkpoint*.mat next to the final file) or cellstr
%     'RelativeThresholds'   [1e-12 1e-6 1e-3 1e-2]
%     'AnalyzeInteraggregate' true
%     'AnalyzePreScatter'    true
%     'MatchPrePost'         true
%     'RunPropertyTests'     true
%     'RenderWorst'          false
%     'MaxRender'            10
%     'SaveResults'          false
%     'OutputDir'            <repoRoot>/results/overlap_audit
%     'Seed'                 20260718 (property-test RNG stream seed)
%     'Verbose'              true
%
%   Example:
%     repoRoot = 'C:\path\to\PFAL';
%     addpath(repoRoot);
%     report = QA.audit_aggregate_overlap(repoRoot, 'AnalyzeInteraggregate', true);
%
%   When repoRoot is omitted or empty it is inferred as the parent folder of
%   the directory containing this file.

% ----------------------------------------------------------------------- %

if nargin < 1 || isempty(repoRoot)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
end
repoRoot = char(repoRoot);

opts = parse_options(repoRoot, varargin{:});

t_all = tic;
report = struct();
report.metadata = struct( ...
    'tool', mfilename('fullpath'), ...
    'repoRoot', repoRoot, ...
    'generated_at', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
    'matlab_version', version);
report.options = opts;

logv(opts, '=== PFAL aggregate overlap / index-integrity audit ===\n');

% ------------------------------------------------------------------- %
% 0. Validate the independent overlap detector
% ------------------------------------------------------------------- %
bp = validate_broadphase(opts);
report.propertyTests.broadphaseValidation = bp;
if ~bp.passed
    error('audit:broadphase', ...
        'Voxel-hash broad phase failed brute-force validation; aborting.');
end
logv(opts, 'Broad-phase self-test: PASSED (%d synthetic + %d real cases)\n', ...
    bp.n_synthetic, bp.n_real);

% ------------------------------------------------------------------- %
% 1. Enumerate stages and load populations one at a time
% ------------------------------------------------------------------- %
[stages, sources, domSize, runMeta] = build_stage_list(opts);
report.sources = sources;

nStg = numel(stages);
stageStats = cell(nStg, 1);
pairRows = cell(nStg, 1);
integRows = cell(nStg, 1);
interRows = cell(nStg, 1);
interStats = cell(nStg, 1);
idmaps = cell(nStg, 1);        % Per-stage [sorted ID, aggregate index] map.
worstPayload = cell(nStg, 1);  % Pair and geometry records used for rendering.

postScatterPP = [];            % Retained for matching, seed mapping, and tests.
seedMap = [];                  % ID-to-seed map; 0 denotes an ambiguous ID.
seedMapMaxId = 0;
collidingLabels = [];          % Column-6 labels shared by multiple seed cells.

for s = 1:nStg
    st = stages(s);
    logv(opts, '\n--- Stage %d/%d: %s (%s) ---\n', s, nStg, st.popId, st.label);
    tS = tic;

    [pp, nStored, nHybStored, extras] = load_stage(st);

    % ---- Internal overlap analysis ---------------------------------- %
    [stat, pairs, worst] = analyze_internal(pp, st, opts);

    % ---- Index-integrity analysis ----------------------------------- %
    integ = integrity_check(pp, nStored, nHybStored, extras, st, opts);

    % ---- Seed map derived from the post-scatter population ---------- %
    if strcmp(st.popId, 'post_scatter')
        postScatterPP = pp;
        [seedMap, seedMapMaxId, collidingLabels] = build_seed_map(pp);
    end

    % ---- Seed-aware overlap and hybridity refinements --------------- %
    if ~isempty(seedMap) && ~strcmp(st.popId, 'pre_scatter')
        [pairs, integ] = refine_with_seeds(pairs, integ, pp, nHybStored, ...
            seedMap, seedMapMaxId, collidingLabels);
    end

    % ---- Interaggregate overlap for LD2 stages ---------------------- %
    if opts.AnalyzeInteraggregate && st.isLD2
        [iStat, iPairs] = analyze_interagg(pp, st, domSize, opts);
        interStats{s} = iStat;
        interRows{s} = iPairs;
    end

    % ---- Primary-ID-to-aggregate map for temporal tracing ----------- %
    if ~strcmp(st.popId, 'pre_scatter')
        idmaps{s} = make_idmap(pp);
    end

    stat.elapsed_s = toc(tS);
    stageStats{s} = stat;
    pairRows{s} = pairs;
    integRows{s} = struct2table(integ, 'AsArray', true);
    worstPayload{s} = worst;

    logv(opts, ['  n_agg=%d  n_pp=%d  pairs(rp>=%g)=%d  maxRelPen=%.3g  ', ...
        '[%.1f s]\n'], stat.nAgg, stat.nPP, opts.RelativeThresholds(1), ...
        stat.pairsAtThr(1), stat.maxRelPen, stat.elapsed_s);

    % Retain the populations required by provenance and rigidity checks.
    if strcmp(st.popId, 'parsdata1')
        parsdata1PP = pp; %#ok<NASGU> % Used by the rigidity check below.
        rigid = post_to_initial_check(postScatterPP, pp, opts);
        report.scatterToLD2InitialCheck = rigid;
    end
    clear pp
end

report.stageSummary = summarize_stages(stageStats, opts);
report.overlapPairs = cat_tables(pairRows);
report.indexIntegrity = cat_tables(integRows);
report.interaggregateSummary = struct( ...
    'perStage', cat_tables(interStats), ...
    'pairs', cat_tables(interRows));

% ------------------------------------------------------------------- %
% 2. Pre-scatter -> post-scatter provenance matching
% ------------------------------------------------------------------- %
if opts.MatchPrePost && ~isempty(postScatterPP)
    report.prePostScatterMatching = match_pre_post(opts, postScatterPP);
else
    report.prePostScatterMatching = struct('performed', false);
end

% ------------------------------------------------------------------- %
% 3. Temporal localization of overlap pairs
% ------------------------------------------------------------------- %
report.checkpointTimeline = build_timeline(stages, stageStats, runMeta, opts);
report.temporalTrace = trace_pairs(stages, pairRows, idmaps, ...
    seedMap, seedMapMaxId, opts);

% ------------------------------------------------------------------- %
% 4. Deterministic collision-function property tests
% ------------------------------------------------------------------- %
if opts.RunPropertyTests
    addpath(repoRoot); % Make PFAL packages available to property-test calls.
    rs = RandStream('twister', 'Seed', opts.Seed);
    report.propertyTests.testA_connect_invariant = test_connect_invariant(postScatterPP, rs, opts);
    report.propertyTests.testB_duplicate_ids = test_duplicate_ids(postScatterPP, opts);
    report.propertyTests.testC_chained_grow = test_chained_grow(opts);
    report.propertyTests.testD_periodic = test_periodic(domSize, opts);
end

% ------------------------------------------------------------------- %
% 5. Worst cases and optional rendering / saving
% ------------------------------------------------------------------- %
report.worstCases = collect_worst(worstPayload, opts);
figs = [];
if opts.RenderWorst
    figs = render_worst(report.worstCases, opts);
end

report.limitations = limitations_text(runMeta, report);

if opts.SaveResults
    report.outputDir = save_outputs(report, figs, opts);
else
    report.outputDir = '';
end

% ------------------------------------------------------------------- %
% 6. Console summary
% ------------------------------------------------------------------- %
print_summary(report, opts);
logv(opts, '\nTotal audit time: %.1f s\n', toc(t_all));

end

% ======================================================================= %
% Configuration and progress reporting
% ======================================================================= %

function opts = parse_options(repoRoot, varargin)
%PARSE_OPTIONS Parse file locations, analysis controls, and output options.
%   Defaults reproduce the archived scatter-to-LD2 run used by this audit.
p = inputParser;
p.FunctionName = 'audit_aggregate_overlap';
dataDir = fullfile(repoRoot, 'data', 'main_scatter');
resDir = fullfile(repoRoot, 'results', 'main_ld2_from_main_scatter', ...
    'LD2__scaled_aggs_for_LD2_from_main_scatter__2026-05-16');
addParameter(p, 'PreScatterFile', fullfile(dataDir, 'LD1__gamma_1_35__merge.mat'));
addParameter(p, 'PostScatterFile', fullfile(dataDir, 'scaled_aggs_for_LD2_from_main_scatter.mat'));
addParameter(p, 'LD2FinalFile', fullfile(resDir, 'LD2__from_main_scatter__final__2026-06-28_20-05-29.mat'));
addParameter(p, 'CheckpointFiles', 'auto');
addParameter(p, 'RelativeThresholds', [1e-12, 1e-6, 1e-3, 1e-2], @(x) isnumeric(x) && all(x > 0));
addParameter(p, 'AnalyzeInteraggregate', true, @islogical);
addParameter(p, 'AnalyzePreScatter', true, @islogical);
addParameter(p, 'MatchPrePost', true, @islogical);
addParameter(p, 'RunPropertyTests', true, @islogical);
addParameter(p, 'RenderWorst', false, @islogical);
addParameter(p, 'MaxRender', 10, @(x) isnumeric(x) && x >= 1);
addParameter(p, 'SaveResults', false, @islogical);
addParameter(p, 'OutputDir', fullfile(repoRoot, 'results', 'overlap_audit'));
addParameter(p, 'Seed', 20260718);
addParameter(p, 'Verbose', true, @islogical);
parse(p, varargin{:});
opts = p.Results;
opts.RelativeThresholds = sort(opts.RelativeThresholds(:)');
opts.repoRoot = repoRoot;
end

function logv(opts, fmt, varargin)
%LOGV Write formatted progress output when verbose reporting is enabled.
if opts.Verbose
    fprintf(fmt, varargin{:});
end
end

% ======================================================================= %
% Stage enumeration and loading
% ======================================================================= %

function [stages, sources, domSize, runMeta] = build_stage_list(opts)
%BUILD_STAGE_LIST Resolve input files and arrange saved populations by time.
%   The returned metadata records incomplete ensemble-history rows and time
%   resets so that gaps introduced by a resumed run remain visible.

% Read run metadata without loading the full final population.
fin = opts.LD2FinalFile;
must_exist(fin);
M = load(fin, 'params_domain', 'n0_agg', 'k', 'ind_dat', 'r_n_agg', 'ensdata');
domSize = M.params_domain.Value(2:4);
domSize = domSize(:)';
runMeta = struct('n0_agg', M.n0_agg, 'k_final', M.k, 'ind_dat_final', M.ind_dat, ...
    'r_n_agg', M.r_n_agg, 'domSize', domSize);
% Derive parsdata save iterations from the trimmed ensemble history. Zero
% rows represent iterations omitted during interruption or resume and must
% be excluded from threshold matching. Associated time resets are retained
% as integrity findings.
nAggHist = M.ensdata.n_agg(:);
zeroRows = find(nAggHist == 0);
kSave = nan(1, numel(M.r_n_agg));
for q = 1:numel(M.r_n_agg)
    kk = find(nAggHist(2:end) <= M.r_n_agg(q) * nAggHist(1) & ...
        nAggHist(2:end) > 0, 1) + 1;
    if ~isempty(kk); kSave(q) = kk; end
end
runMeta.parsdata_k = kSave;
runMeta.nAggHist = nAggHist;
runMeta.ensdataZeroRows = zeroRows;
tHist = M.ensdata.t(:);
runMeta.ensdataTimeResets = find(diff(tHist) < 0) + 1; % Cumulative-time resets.
if ~isempty(zeroRows)
    fprintf(['NOTE: ensdata.n_agg has %d unrecorded (zero) row(s) at k=[%s]; ', ...
        'cumulative time resets at k=[%s] -- evidence of mid-run resume.\n'], ...
        numel(zeroRows), num2str(zeroRows(:)'), ...
        num2str(runMeta.ensdataTimeResets(:)'));
end

% Discover checkpoint files unless an explicit list was supplied.
if ischar(opts.CheckpointFiles) || isstring(opts.CheckpointFiles)
    if strcmpi(char(opts.CheckpointFiles), 'auto')
        d = dir(fullfile(fileparts(fin), '*checkpoint*.mat'));
        cpFiles = sort(fullfile({d.folder}, {d.name}));
    else
        cpFiles = {char(opts.CheckpointFiles)};
    end
else
    cpFiles = cellfun(@char, opts.CheckpointFiles, 'UniformOutput', false);
end

stages = struct('popId', {}, 'label', {}, 'file', {}, 'kind', {}, ...
    'slot', {}, 'isLD2', {}, 'kIter', {});
add = @(st, popId, label, file, kind, slot, isLD2, kIter) ...
    [st, struct('popId', popId, 'label', label, 'file', file, 'kind', kind, ...
    'slot', slot, 'isLD2', isLD2, 'kIter', kIter)];

if opts.AnalyzePreScatter
    must_exist(opts.PreScatterFile);
    stages = add(stages, 'pre_scatter', 'candidate raw LD1 library pp0', ...
        opts.PreScatterFile, 'pp0', 0, false, NaN);
end
must_exist(opts.PostScatterFile);
stages = add(stages, 'post_scatter', 'scatter output pars_out (LD2 input)', ...
    opts.PostScatterFile, 'pars_out', 0, false, NaN);
for q = 1:5
    stages = add(stages, sprintf('parsdata%d', q), ...
        sprintf('LD2 saved population %d', q), fin, 'parsdata', q, true, kSave(q));
end
for q = 1:numel(cpFiles)
    must_exist(cpFiles{q});
    C = load(cpFiles{q}, 'k', 'ind_dat');
    stages = add(stages, sprintf('checkpoint%d', q), ...
        sprintf('LD2 checkpoint (k=%d, ind_dat=%d)', C.k, C.ind_dat), ...
        cpFiles{q}, 'pars_LD2', q, true, C.k);
end
stages = add(stages, 'final', 'final pars_LD2', fin, 'pars_LD2', 0, true, M.k);

% Place source populations first and order LD2 populations by iteration.
kOrd = [stages.kIter];
kOrd(isnan(kOrd)) = -inf;
base = zeros(1, numel(stages));
base(strcmp({stages.popId}, 'pre_scatter')) = -2;
base(strcmp({stages.popId}, 'post_scatter')) = -1;
[~, ord] = sortrows([base(:), kOrd(:)]);
stages = stages(ord);

files = unique({stages.file});
sources = struct('files', {files});
for q = 1:numel(files)
    d = dir(files{q});
    sources.bytes(q) = d.bytes;
end
end

function must_exist(f)
%MUST_EXIST Raise a descriptive error when a required input is unavailable.
if ~isfile(f)
    error('audit:missingFile', 'Required data file not found: %s', f);
end
end

function [pp, nStored, nHybStored, extras] = load_stage(st)
%LOAD_STAGE Load one population and the stored fields required by the audit.
%   Loading one stage at a time limits peak memory use for large MAT files.
nStored = [];
nHybStored = [];
extras = struct();
switch st.kind
    case 'pp0'
        S = load(st.file, 'pp0');
        pp = S.pp0;
    case 'pars_out'
        S = load(st.file, 'pars_out');
        pp = S.pars_out.pp;
        if isfield(S.pars_out, 'n'); nStored = S.pars_out.n; end
        fn = intersect(fieldnames(S.pars_out), ...
            {'n', 'dpp_g', 'da', 'dv', 'dg', 'dmax', 'dpp'});
        for q = 1:numel(fn)
            extras.(fn{q}) = S.pars_out.(fn{q});
        end
    case 'parsdata'
        S = load(st.file, 'parsdata');
        P = S.parsdata(st.slot);
        pp = P.pp;
        nStored = P.npp;
        if isfield(P, 'n_hyb'); nHybStored = P.n_hyb; end
        fn = intersect(fieldnames(P), {'dpp', 'sigmapp', 'da', 'dm', 'dg'});
        for q = 1:numel(fn)
            extras.(fn{q}) = P.(fn{q});
        end
    case 'pars_LD2'
        S = load(st.file, 'pars_LD2');
        pp = S.pars_LD2.pp;
        nStored = S.pars_LD2.n;
        if isfield(S.pars_LD2, 'n_hyb'); nHybStored = S.pars_LD2.n_hyb; end
        fn = intersect(fieldnames(S.pars_LD2), ...
            {'n', 'r', 'v', 'm', 'dmax', 'dpp_g', 'da', 'dv', 'dg', 'dm', 'dpp'});
        for q = 1:numel(fn)
            extras.(fn{q}) = S.pars_LD2.(fn{q});
        end
end
pp = pp(:);
end

% ======================================================================= %
% Independent overlap detection (voxel hash + brute force)
% ======================================================================= %

function cand = candidate_pairs_self(xyz, d)
%CANDIDATE_PAIRS_SELF Find candidate pairs within one particle population.
%   The broad phase returns each index pair [i,j], i < j, whose center
%   distance is less than max(d). This is a superset of all overlapping
%   pairs because no pair has a contact distance greater than max(d).
n = size(xyz, 1);
if n < 2
    cand = zeros(0, 2);
elseif n <= 350
    cand = brute_candidates_self(xyz, max(d));
else
    cand = grid_candidates_self(xyz, max(d));
end
end

function cand = brute_candidates_self(xyz, cutoff)
%BRUTE_CANDIDATES_SELF Enumerate self-pairs by blocked distance evaluation.
%   Blocking bounds temporary matrix size while preserving an exact search.
n = size(xyz, 1);
cand = cell(0, 1);
blk = 512;
c2 = cutoff^2;
for i0 = 1:blk:n
    i1 = min(i0 + blk - 1, n);
    ii = (i0:i1)';
    dx = xyz(ii, 1) - xyz(:, 1)';
    dy = xyz(ii, 2) - xyz(:, 2)';
    dz = xyz(ii, 3) - xyz(:, 3)';
    d2 = dx.^2 + dy.^2 + dz.^2;
    mask = d2 < c2 & (ii < (1:n)); % Retain each unordered pair once.
    [ra, ca] = find(mask);
    cand{end+1, 1} = [ii(ra), ca]; %#ok<AGROW>
end
cand = vertcat(cand{:});
if isempty(cand); cand = zeros(0, 2); end
end

function cand = grid_candidates_self(xyz, cutoff)
%GRID_CANDIDATES_SELF Enumerate self-pairs with an exact voxel hash.
%   A voxel width equal to the cutoff guarantees that every qualifying pair
%   occupies either the same voxel or one of its 26 adjacent voxels.
n = size(xyz, 1);
g = floor(xyz ./ cutoff);
g = g - min(g, [], 1) + 1;               % Shift grid coordinates to positive indices.
dims = max(g, [], 1);
M1 = dims(1) + 3; M2 = dims(2) + 3;      % Define an alias-free mixed-radix key.
key = g(:, 1) + M1 .* (g(:, 2) + M2 .* g(:, 3));
[uk, ~, ic] = unique(key);
[sic, ord] = sort(ic);
b0 = [1; find(diff(sic)) + 1];           % Locate each bucket start in ord.
b1 = [b0(2:end) - 1; n];
% Use 13 half-space offsets so that each adjacent bucket pair is visited once.
offs = [];
for oz = -1:1
    for oy = -1:1
        for ox = -1:1
            if (oz > 0) || (oz == 0 && oy > 0) || (oz == 0 && oy == 0 && ox > 0)
                offs(end+1, :) = [ox, oy, oz]; %#ok<AGROW>
            end
        end
    end
end
okey = offs(:, 1) + M1 .* (offs(:, 2) + M2 .* offs(:, 3));
cand = cell(0, 1);
% Enumerate pairs whose particles occupy the same voxel.
for b = 1:numel(uk)
    m = ord(b0(b):b1(b));
    if numel(m) > 1
        pr = nchoosek(sort(m), 2);
        cand{end+1, 1} = pr; %#ok<AGROW>
    end
end
% Enumerate pairs in adjacent voxels.
for o = 1:numel(okey)
    [tf, loc] = ismember(uk + okey(o), uk);
    src = find(tf);
    for q = 1:numel(src)
        a = ord(b0(src(q)):b1(src(q)));
        b = ord(b0(loc(src(q))):b1(loc(src(q))));
        [A, B] = ndgrid(a, b);
        pr = [A(:), B(:)];
        swap = pr(:, 1) > pr(:, 2);
        pr(swap, :) = pr(swap, [2, 1]);
        cand{end+1, 1} = pr; %#ok<AGROW>
    end
end
cand = vertcat(cand{:});
if isempty(cand)
    cand = zeros(0, 2);
else
    cand = unique(cand, 'rows'); % Remove any duplicate candidate indices.
end
end

function cand = candidate_pairs_cross(xyz1, xyz2, cutoff)
%CANDIDATE_PAIRS_CROSS Find candidate pairs between two populations.
%   Returns [i,j] when the distance between xyz1(i,:) and xyz2(j,:) is less
%   than the supplied cutoff. Direct blocking is used for smaller products;
%   larger products use the same exact voxel-neighborhood guarantee.
n1 = size(xyz1, 1);
n2 = size(xyz2, 1);
if n1 * n2 <= 250000
    cand = cell(0, 1);
    c2 = cutoff^2;
    blk = 512;
    for i0 = 1:blk:n1
        i1 = min(i0 + blk - 1, n1);
        ii = (i0:i1)';
        d2 = (xyz1(ii, 1) - xyz2(:, 1)').^2 + (xyz1(ii, 2) - xyz2(:, 2)').^2 ...
            + (xyz1(ii, 3) - xyz2(:, 3)').^2;
        [ra, ca] = find(d2 < c2);
        cand{end+1, 1} = [ii(ra), ca]; %#ok<AGROW>
    end
    cand = vertcat(cand{:});
    if isempty(cand); cand = zeros(0, 2); end
    return
end
% Hash xyz2 and query xyz1 over the complete 27-voxel neighborhood.
orig = min([xyz1; xyz2], [], 1);
g2 = floor((xyz2 - orig) ./ cutoff) + 1;
g1 = floor((xyz1 - orig) ./ cutoff) + 1;
dims = max([g1; g2], [], 1);
M1 = dims(1) + 3; M2 = dims(2) + 3;
key2 = g2(:, 1) + M1 .* (g2(:, 2) + M2 .* g2(:, 3));
[uk, ~, ic] = unique(key2);
[sic, ord] = sort(ic);
b0 = [1; find(diff(sic)) + 1];
b1 = [b0(2:end) - 1; n2];
key1 = g1(:, 1) + M1 .* (g1(:, 2) + M2 .* g1(:, 3));
cand = cell(0, 1);
for ox = -1:1
    for oy = -1:1
        for oz = -1:1
            ok = ox + M1 .* (oy + M2 .* oz);
            [tf, loc] = ismember(key1 + ok, uk);
            src = find(tf);
            for q = 1:numel(src)
                b = ord(b0(loc(src(q))):b1(loc(src(q))));
                cand{end+1, 1} = [repmat(src(q), numel(b), 1), b(:)]; %#ok<AGROW>
            end
        end
    end
end
cand = vertcat(cand{:});
if isempty(cand); cand = zeros(0, 2); end
end

function res = exact_pairs(xyz1, d1, xyz2, d2, cand)
%EXACT_PAIRS Apply the geometric overlap test to broad-phase candidates.
%   Columns are [row1,row2,center distance,contact distance,clearance]. A
%   negative clearance denotes overlap. Passing the same arrays in both
%   positions performs a self-population narrow phase.
if isempty(cand)
    res = zeros(0, 5);
    return
end
dv = xyz1(cand(:, 1), :) - xyz2(cand(:, 2), :);
dist = sqrt(sum(dv.^2, 2));
contact = (d1(cand(:, 1)) + d2(cand(:, 2))) ./ 2;
res = [cand, dist, contact, dist - contact]; % Preserve geometry for downstream reporting.
end

function bp = validate_broadphase(opts)
%VALIDATE_BROADPHASE Compare voxel-hash results with brute-force results.
%   Validation covers seeded synthetic clouds at multiple length scales,
%   cross-population searches, and a sample of stored post-scatter
%   aggregates. Any disagreement invalidates the detector.
rs = RandStream('twister', 'Seed', 12345);
bp = struct('passed', true, 'n_synthetic', 0, 'n_real', 0, 'failures', {{}});
for trial = 1:20
    n = 30 + randi(rs, 400);
    scale = 10^(randi(rs, [3, 9]) * -1); % Exercise multiple absolute length scales.
    xyz = rand(rs, n, 3) .* scale;
    d = (0.02 + 0.2 * rand(rs, n, 1)) .* scale .* (1 + 4 * rand(rs, n, 1));
    % Construct the full O(n^2) reference pair set.
    ref = zeros(0, 2);
    for i = 1:n-1
        dd = sqrt(sum((xyz(i+1:end, :) - xyz(i, :)).^2, 2));
        j = find(dd < (d(i) + d(i+1:end)) ./ 2) + i;
        ref = [ref; repmat(i, numel(j), 1), j(:)]; %#ok<AGROW>
    end
    % Force the voxel-hash path and compare exact overlap sets.
    cand = grid_candidates_self(xyz, max(d));
    got = exact_pairs(xyz, d, xyz, d, cand);
    got = sortrows(got(got(:, 5) < 0, 1:2));
    if ~isequal(sortrows(ref), got)
        bp.passed = false;
        bp.failures{end+1} = sprintf('synthetic trial %d', trial);
    end
    bp.n_synthetic = bp.n_synthetic + 1;
    % Validate the cross-population path by splitting the same cloud.
    h = floor(n / 2);
    candX = candidate_pairs_cross(xyz(1:h, :), xyz(h+1:end, :), max(d));
    gotX = exact_pairs(xyz(1:h, :), d(1:h), xyz(h+1:end, :), d(h+1:end), candX);
    gotX = sortrows(gotX(gotX(:, 5) < 0, 1:2));
    refX = ref(ref(:, 1) <= h & ref(:, 2) > h, :);
    refX = sortrows([refX(:, 1), refX(:, 2) - h]);
    if ~isequal(refX, gotX)
        bp.passed = false;
        bp.failures{end+1} = sprintf('synthetic cross trial %d', trial);
    end
end
% Repeat the comparison on stored aggregates small enough for brute force.
try
    S = load(opts.PostScatterFile, 'pars_out');
    ppl = S.pars_out.pp;
    sizes = cellfun(@(x) size(x, 1), ppl);
    pick = find(sizes >= 5 & sizes <= 500);
    pick = pick(round(linspace(1, numel(pick), min(30, numel(pick)))));
    for q = 1:numel(pick)
        P = ppl{pick(q)};
        xyz = P(:, 3:5); d = P(:, 2);
        cand = grid_candidates_self(xyz, max(d));
        got = exact_pairs(xyz, d, xyz, d, cand);
        got = sortrows(got(got(:, 5) < 0, 1:2));
        ref = brute_candidates_self(xyz, max(d));
        refE = exact_pairs(xyz, d, xyz, d, ref);
        refE = sortrows(refE(refE(:, 5) < 0, 1:2));
        if ~isequal(ref_all(refE), ref_all(got))
            bp.passed = false;
            bp.failures{end+1} = sprintf('real aggregate %d', pick(q));
        end
        bp.n_real = bp.n_real + 1;
    end
catch err
    bp.failures{end+1} = ['real-aggregate check skipped: ', err.message];
end
end

function x = ref_all(x)
%REF_ALL Preserve a two-column shape for empty reference pair sets.
if isempty(x); x = zeros(0, 2); end
end

% ======================================================================= %
% Per-stage internal overlap analysis
% ======================================================================= %

function [stat, pairsT, worst] = analyze_internal(pp, st, opts)
%ANALYZE_INTERNAL Measure primary-particle overlap within each aggregate.
%   Pair records retain particle identity, labels, geometry, and penetration
%   depth. Summary counts are evaluated independently at every configured
%   relative-penetration threshold.
thr = opts.RelativeThresholds;
nAgg = numel(pp);
nppv = cellfun(@(x) size(x, 1), pp);
recThr = thr(1);

recs = cell(nAgg, 1);
minClr = inf;
nContact = 0;
nRoundoff = 0;
for i = 1:nAgg
    P = pp{i};
    if size(P, 1) < 2; continue; end
    xyz = P(:, 3:5); d = P(:, 2);
    cand = candidate_pairs_self(xyz, d);
    res = exact_pairs(xyz, d, xyz, d, cand);
    if isempty(res); continue; end
    minClr = min(minClr, min(res(:, 5)));
    relp = -res(:, 5) ./ res(:, 4);
    nContact = nContact + nnz(abs(relp) <= 1e-12);
    nRoundoff = nRoundoff + nnz(relp > 0 & relp <= recThr);
    keep = relp > recThr;
    if any(keep)
        r = res(keep, :);
        rp = relp(keep);
        ids = P(:, 1); lbl = P(:, 6);
        dupIds = ids(count_dups(ids));
        rec = struct();
        rec.aggIdx = repmat(i, size(r, 1), 1);
        rec.row1 = r(:, 1); rec.row2 = r(:, 2);
        rec.id1 = ids(r(:, 1)); rec.id2 = ids(r(:, 2));
        rec.lbl1 = lbl(r(:, 1)); rec.lbl2 = lbl(r(:, 2));
        rec.d1 = d(r(:, 1)); rec.d2 = d(r(:, 2));
        rec.x1 = xyz(r(:, 1), 1); rec.y1 = xyz(r(:, 1), 2); rec.z1 = xyz(r(:, 1), 3);
        rec.x2 = xyz(r(:, 2), 1); rec.y2 = xyz(r(:, 2), 2); rec.z2 = xyz(r(:, 2), 3);
        rec.centerDist = r(:, 3); rec.contactDist = r(:, 4);
        rec.clearance = r(:, 5);
        rec.penetration = -r(:, 5);
        rec.relPen = rp;
        rec.idsEqual = rec.id1 == rec.id2;
        rec.lblsEqual = rec.lbl1 == rec.lbl2;
        rec.id1DupInAgg = ismember(rec.id1, dupIds);
        rec.id2DupInAgg = ismember(rec.id2, dupIds);
        recs{i} = rec;
    end
end
recs = recs(~cellfun(@isempty, recs));
if isempty(recs)
    pairsT = table();
else
    pairsT = struct2table(merge_recs(recs));
    pairsT.stage = repmat({st.popId}, height(pairsT), 1);
    pairsT.kIter = repmat(st.kIter, height(pairsT), 1);
    pairsT.seed1 = nan(height(pairsT), 1);
    pairsT.seed2 = nan(height(pairsT), 1);
    pairsT.lblCollides = false(height(pairsT), 1);
    pairsT = movevars_safe(pairsT, 'stage', 1);
end

stat = struct();
stat.popId = st.popId;
stat.label = st.label;
stat.file = st.file;
stat.kIter = st.kIter;
stat.nAgg = nAgg;
stat.nPP = sum(nppv);
stat.minAggSize = min(nppv);
stat.maxAggSize = max(nppv);
stat.minClearance = minClr;
stat.contactPairs = nContact;
stat.roundoffOnlyPairs = nRoundoff;
nT = numel(thr);
stat.pairsAtThr = zeros(1, nT);
stat.aggsAtThr = zeros(1, nT);
stat.ppAtThr = zeros(1, nT);
if ~isempty(pairsT)
    for q = 1:nT
        m = pairsT.relPen >= thr(q);
        stat.pairsAtThr(q) = nnz(m);
        stat.aggsAtThr(q) = numel(unique(pairsT.aggIdx(m)));
        k1 = [pairsT.aggIdx(m), pairsT.row1(m)];
        k2 = [pairsT.aggIdx(m), pairsT.row2(m)];
        stat.ppAtThr(q) = size(unique([k1; k2], 'rows'), 1);
    end
    stat.maxAbsPen = max(pairsT.penetration);
    stat.maxRelPen = max(pairsT.relPen);
    posrp = pairsT.relPen;
    stat.medianRelPen = median(posrp);
    stat.p95RelPen = local_prctile(posrp, 95);
    stat.p99RelPen = local_prctile(posrp, 99);
    stat.sameLblPairs = nnz(pairsT.lblsEqual);
    stat.crossLblPairs = nnz(~pairsT.lblsEqual);
    stat.dupIdPairs = nnz(pairsT.id1DupInAgg | pairsT.id2DupInAgg | pairsT.idsEqual);
else
    [stat.maxAbsPen, stat.maxRelPen, stat.medianRelPen, ...
        stat.p95RelPen, stat.p99RelPen] = deal(0);
    [stat.sameLblPairs, stat.crossLblPairs, stat.dupIdPairs] = deal(0);
end
stat.dupLblPairs = NaN; % Populated after source-seed labels are available.

% Retain the most severe pairs and their geometry for optional rendering.
worst = struct('stage', {}, 'kIter', {}, 'aggIdx', {}, 'row1', {}, 'row2', {}, ...
    'relPen', {}, 'clearance', {}, 'pp', {});
if ~isempty(pairsT)
    [~, ordw] = sort(pairsT.relPen, 'descend');
    ordw = ordw(1:min(10, numel(ordw)));
    for q = 1:numel(ordw)
        rw = pairsT(ordw(q), :);
        worst(q).stage = st.popId;
        worst(q).kIter = st.kIter;
        worst(q).aggIdx = rw.aggIdx;
        worst(q).row1 = rw.row1;
        worst(q).row2 = rw.row2;
        worst(q).relPen = rw.relPen;
        worst(q).clearance = rw.clearance;
        worst(q).pp = pp{rw.aggIdx};
    end
end
end

function tf = count_dups(v)
%COUNT_DUPS Return a mask for values that occur more than once.
[sv, ~] = sort(v);
dupVals = sv([false; diff(sv) == 0]);
tf = ismember(v, dupVals);
end

function m = merge_recs(recs)
%MERGE_RECS Concatenate fields from a cell array of scalar record structs.
fn = fieldnames(recs{1});
m = struct();
for q = 1:numel(fn)
    parts = cellfun(@(r) r.(fn{q}), recs, 'UniformOutput', false);
    m.(fn{q}) = vertcat(parts{:});
end
end

function T = movevars_safe(T, name, pos)
%MOVEVARS_SAFE Reorder a table variable when movevars is available.
%   Table contents remain valid on MATLAB releases that lack movevars.
try
    T = movevars(T, name, 'Before', pos);
catch
end
end

function v = local_prctile(x, p)
%LOCAL_PRCTILE Compute a linearly interpolated percentile without toolboxes.
x = sort(x(:));
n = numel(x);
if n == 0; v = NaN; return; end
if n == 1; v = x; return; end
q = (p / 100) * (n - 1) + 1;
lo = floor(q); hi = ceil(q);
v = x(lo) + (q - lo) * (x(hi) - x(lo));
end

% ======================================================================= %
% Index integrity
% ======================================================================= %

function integ = integrity_check(pp, nStored, nHybStored, extras, st, opts) %#ok<INUSD>
%INTEGRITY_CHECK Compare particle arrays with stored identifiers and metadata.
%   The checks cover primary IDs, subaggregate labels, aggregate-level array
%   lengths, stored centers of mass, bounding diameters, and hybrid counts.
nAgg = numel(pp);
nppv = cellfun(@(x) size(x, 1), pp);

ids = cellfun(@(x) x(:, 1), pp, 'UniformOutput', false);
lbl = cellfun(@(x) x(:, 6), pp, 'UniformOutput', false);
allIds = vertcat(ids{:});
allLbl = vertcat(lbl{:});

integ = struct();
integ.stage = {st.popId};
integ.nAgg = nAgg;
integ.nPP = sum(nppv);

% Compare stored primary-particle counts with actual matrix row counts.
if isempty(nStored)
    integ.nStoredMismatch = NaN;
else
    integ.nStoredMismatch = nnz(nStored(:) ~= nppv(:));
end

% Audit primary-particle identifiers stored in column 1.
integ.id_nonfinite = nnz(~isfinite(allIds));
integ.id_noninteger = nnz(isfinite(allIds) & allIds ~= round(allIds));
integ.id_nonpositive = nnz(allIds <= 0);
dupWithin = 0;
aggsWithDupIds = 0;
for i = 1:nAgg
    ndup = nppv(i) - numel(unique(ids{i}));
    dupWithin = dupWithin + ndup;
    aggsWithDupIds = aggsWithDupIds + (ndup > 0);
end
integ.id_dupWithinAgg = dupWithin;
integ.id_aggsWithDupIds = aggsWithDupIds;
integ.id_dupAcrossCells = numel(allIds) - numel(unique(allIds));

% Audit subaggregate labels stored in column 6.
integ.lbl_nonfinite = nnz(~isfinite(allLbl));
integ.lbl_noninteger = nnz(isfinite(allLbl) & allLbl ~= round(allLbl));
integ.lbl_nonpositive = nnz(allLbl <= 0);
firstLbl = cellfun(@(x) x(1), lbl);
uLblPerAgg = cellfun(@(x) numel(unique(x)), lbl);
% Count labels shared by multiple aggregate cells. This metric is meaningful
% only while each cell still represents a single source seed.
if all(uLblPerAgg == 1)
    integ.lbl_cellsSharingLabel = nnz(count_dups(firstLbl));
    integ.lbl_distinctSharedLabels = numel(unique(firstLbl(count_dups(firstLbl))));
else
    integ.lbl_cellsSharingLabel = NaN;
    integ.lbl_distinctSharedLabels = NaN;
end

% Compare stored hybrid counts with the number of distinct column-6 labels.
if isempty(nHybStored)
    integ.nHyb_mismatch = NaN;
    integ.nHyb_max = NaN;
else
    integ.nHyb_mismatch = nnz(nHybStored(:) ~= uLblPerAgg(:));
    integ.nHyb_max = max(nHybStored);
end
integ.nHyb_trueUndercount = NaN; % Populated after source-seed mapping.

% Check that every aggregate-level array has one row per aggregate cell.
fn = fieldnames(extras);
mis = {};
for q = 1:numel(fn)
    v = extras.(fn{q});
    if size(v, 1) ~= nAgg
        mis{end+1} = sprintf('%s:%d', fn{q}, size(v, 1)); %#ok<AGROW>
    end
end
integ.arrayLenMismatch = {strjoin(mis, ',')};

% Validate stored centers of mass and bounding diameters when available.
integ.max_r_vs_com = NaN;
integ.max_dmax_relErr = NaN;
if isfield(extras, 'r') && size(extras.r, 1) == nAgg
    dev = zeros(nAgg, 1);
    for i = 1:nAgg
        w = pp{i}(:, 2).^3;
        com = sum(w .* pp{i}(:, 3:5), 1) ./ sum(w);
        dev(i) = norm(extras.r(i, :) - com);
    end
    integ.max_r_vs_com = max(dev);
end
if isfield(extras, 'dmax') && size(extras.dmax, 1) == nAgg
    rel = zeros(nAgg, 1);
    for i = 1:nAgg
        w = pp{i}(:, 2).^3;
        com = sum(w .* pp{i}(:, 3:5), 1) ./ sum(w);
        dm = max(2 .* sqrt(sum((pp{i}(:, 3:5) - com).^2, 2)) + pp{i}(:, 2));
        rel(i) = abs(extras.dmax(i) - dm) / dm;
    end
    integ.max_dmax_relErr = max(rel);
end
end

function [seedMap, maxId, collidingLabels] = build_seed_map(pp)
%BUILD_SEED_MAP Map unambiguous primary IDs to post-scatter source cells.
%   IDs appearing in more than one source cell receive a zero mapping and
%   are excluded from provenance claims that require unique identity.
ids = cellfun(@(x) x(:, 1), pp, 'UniformOutput', false);
allIds = vertcat(ids{:});
maxId = max(allIds);
seedMap = zeros(maxId, 1, 'uint32');
amb = false(maxId, 1);
for i = 1:numel(pp)
    v = ids{i};
    ok = v >= 1 & v == round(v);
    v = v(ok);
    amb(v(seedMap(v) ~= 0)) = true;
    seedMap(v) = i;
end
seedMap(amb) = 0; % Exclude ambiguous IDs from source-seed tracing.
lblFirst = cellfun(@(x) x(1, 6), pp);
collidingLabels = unique(lblFirst(count_dups(lblFirst)));
end

function [pairsT, integ] = refine_with_seeds(pairsT, integ, pp, nHybStored, ...
    seedMap, maxId, collidingLabels)
%REFINE_WITH_SEEDS Add provenance fields and seed-based hybridity checks.
%   Pair-level source identities use only unambiguous primary IDs. Aggregate
%   hybrid counts compare distinct source seeds with distinct stored labels.

% Add source-seed identity and label-collision status to each overlap pair.
if ~isempty(pairsT)
    s1 = zeros(height(pairsT), 1); s2 = zeros(height(pairsT), 1);
    ok1 = pairsT.id1 >= 1 & pairsT.id1 <= maxId & pairsT.id1 == round(pairsT.id1);
    ok2 = pairsT.id2 >= 1 & pairsT.id2 <= maxId & pairsT.id2 == round(pairsT.id2);
    s1(ok1) = double(seedMap(pairsT.id1(ok1)));
    s2(ok2) = double(seedMap(pairsT.id2(ok2)));
    pairsT.seed1 = s1;
    pairsT.seed2 = s2;
    pairsT.lblCollides = ismember(pairsT.lbl1, collidingLabels) | ...
        ismember(pairsT.lbl2, collidingLabels);
end
% Compare source-seed hybridity with the label-based stored hybrid count.
if ~isempty(nHybStored)
    nAgg = numel(pp);
    under = 0;
    for i = 1:nAgg
        v = pp{i}(:, 1);
        ok = v >= 1 & v <= maxId & v == round(v);
        sd = double(seedMap(v(ok)));
        sd = sd(sd > 0);
        trueHyb = numel(unique(sd));
        lblHyb = numel(unique(pp{i}(:, 6)));
        under = under + max(0, trueHyb - lblHyb);
    end
    integ.nHyb_trueUndercount = under;
end
end

function map = make_idmap(pp)
%MAKE_IDMAP Build a sorted primary-ID-to-aggregate lookup for one stage.
%   Duplicate primary IDs remain duplicated in the lookup. Consequently,
%   temporal traces keyed only by ID pairs are not uniquely attributable for
%   the ambiguous subset; geometric overlap counts do not use this lookup.
ids = cellfun(@(x) x(:, 1), pp, 'UniformOutput', false);
agg = cell(numel(pp), 1);
for i = 1:numel(pp)
    agg{i} = repmat(i, numel(ids{i}), 1);
end
allIds = vertcat(ids{:});
allAgg = vertcat(agg{:});
[map.id, ord] = sort(allIds);
map.agg = allAgg(ord);
end

% ======================================================================= %
% Interaggregate overlap
% ======================================================================= %

function [statT, pairsT] = analyze_interagg(pp, st, domSize, opts)
%ANALYZE_INTERAGG Measure overlap between distinct aggregate cells.
%   The ordinary mode uses stored coordinates directly. The periodic mode
%   also evaluates relevant periodic images under the minimum-image domain.
thr = opts.RelativeThresholds;
recThr = thr(1);
nAgg = numel(pp);
com = zeros(nAgg, 3);
R = zeros(nAgg, 1);
maxd = zeros(nAgg, 1);
for i = 1:nAgg
    w = pp{i}(:, 2).^3;
    com(i, :) = sum(w .* pp{i}(:, 3:5), 1) ./ sum(w);
    R(i) = max(sqrt(sum((pp{i}(:, 3:5) - com(i, :)).^2, 2)) + pp{i}(:, 2) ./ 2);
    maxd(i) = max(pp{i}(:, 2));
end

modes = {'ordinary', 'periodic'};
rows = cell(2, 1);
statRows = cell(2, 1);
for md = 1:2
    periodic = strcmp(modes{md}, 'periodic');
    % Use aggregate bounding spheres to form a chunked candidate-pair set.
    candAgg = zeros(0, 2);
    blk = 256;
    for i0 = 1:blk:nAgg
        i1 = min(i0 + blk - 1, nAgg);
        ii = (i0:i1)';
        dx = com(ii, 1) - com(:, 1)';
        dy = com(ii, 2) - com(:, 2)';
        dz = com(ii, 3) - com(:, 3)';
        if periodic
            dx = dx - domSize(1) .* round(dx ./ domSize(1));
            dy = dy - domSize(2) .* round(dy ./ domSize(2));
            dz = dz - domSize(3) .* round(dz ./ domSize(3));
        end
        dc = sqrt(dx.^2 + dy.^2 + dz.^2);
        lim = R(ii) + R' + max(maxd(ii), maxd'); % Include the primary-pair cutoff.
        mask = dc < lim & (ii < (1:nAgg));
        [ra, ca] = find(mask);
        candAgg = [candAgg; ii(ra), ca]; %#ok<AGROW>
    end
    % Apply the exact primary-particle test to each aggregate candidate.
    rec = cell(size(candAgg, 1), 1);
    for q = 1:size(candAgg, 1)
        a = candAgg(q, 1); b = candAgg(q, 2);
        P1 = pp{a}; P2 = pp{b};
        cutoff = (maxd(a) + maxd(b)) / 2;
        if periodic
            imgs = image_shifts(com(a, :), com(b, :), R(a) + R(b) + cutoff, domSize);
        else
            imgs = [0, 0, 0];
        end
        for gsh = 1:size(imgs, 1)
            sh = imgs(gsh, :);
            xyz2 = P2(:, 3:5) + sh;
            cand = candidate_pairs_cross(P1(:, 3:5), xyz2, cutoff);
            res = exact_pairs(P1(:, 3:5), P1(:, 2), xyz2, P2(:, 2), cand);
            if isempty(res); continue; end
            relp = -res(:, 5) ./ res(:, 4);
            keep = relp > recThr;
            if ~any(keep); continue; end
            r = res(keep, :);
            c = struct();
            c.agg1 = repmat(a, size(r, 1), 1);
            c.agg2 = repmat(b, size(r, 1), 1);
            c.row1 = r(:, 1); c.row2 = r(:, 2);
            c.id1 = P1(r(:, 1), 1); c.id2 = P2(r(:, 2), 1);
            c.lbl1 = P1(r(:, 1), 6); c.lbl2 = P2(r(:, 2), 6);
            c.centerDist = r(:, 3); c.contactDist = r(:, 4);
            c.clearance = r(:, 5); c.penetration = -r(:, 5);
            c.relPen = relp(keep);
            c.shiftX = repmat(sh(1), size(r, 1), 1);
            c.shiftY = repmat(sh(2), size(r, 1), 1);
            c.shiftZ = repmat(sh(3), size(r, 1), 1);
            rec{q} = [rec{q}; struct2table(c)]; %#ok<AGROW>
        end
    end
    rec = rec(~cellfun(@isempty, rec));
    if isempty(rec)
        T = table();
    else
        T = vertcat(rec{:});
        T.stage = repmat({st.popId}, height(T), 1);
        T.mode = repmat(modes(md), height(T), 1);
        T = movevars_safe(T, 'stage', 1);
    end
    rows{md} = T;

    srow = struct();
    srow.stage = {st.popId};
    srow.mode = modes(md);
    srow.kIter = st.kIter;
    srow.candidateAggPairs = size(candAgg, 1);
    if isempty(T)
        srow.overlappingAggPairs = 0;
        srow.primaryPairs = 0;
        srow.maxRelPen = 0;
        srow.materialPairs = 0;
    else
        srow.overlappingAggPairs = size(unique([T.agg1, T.agg2], 'rows'), 1);
        srow.primaryPairs = height(T);
        srow.maxRelPen = max(T.relPen);
        srow.materialPairs = nnz(T.relPen >= 1e-3);
    end
    statRows{md} = struct2table(srow);
end
pairsT = cat_tables(rows);
statT = vertcat(statRows{:});
end

function imgs = image_shifts(com1, com2, reach, domSize)
%IMAGE_SHIFTS Return periodic images that may fall within the search reach.
%   The nearest image is always considered. Adjacent images are included
%   when their aggregate bounding regions can also intersect.
base = -domSize .* round((com2 - com1) ./ domSize);
imgs = zeros(0, 3);
for ox = -1:1
    for oy = -1:1
        for oz = -1:1
            sh = base + [ox, oy, oz] .* domSize;
            if norm(com2 + sh - com1) <= reach
                imgs(end+1, :) = sh; %#ok<AGROW>
            end
        end
    end
end
if isempty(imgs); imgs = base; end
end

% ======================================================================= %
% Pre-scatter -> post-scatter matching
% ======================================================================= %

function M = match_pre_post(opts, postPP)
%MATCH_PRE_POST Match post-scatter aggregates to their pre-scatter sources.
%   Candidate matches use a compact ID-multiset fingerprint and are then
%   verified by exact sorted-ID equality. Matched geometry is tested for a
%   uniform diameter scale, translation, and preservation of pair distances.
logv(opts, '\n--- Pre->post scatter provenance matching ---\n');
S = load(opts.PreScatterFile, 'pp0');
prePP = S.pp0(:);
clear S
nPre = numel(prePP);
nPost = numel(postPP);

% Index the pre-scatter library by an ID-multiset fingerprint.
keyOf = @(idv) sprintf('%d|%.0f|%.0f|%.0f', numel(idv), min(idv), max(idv), sum(idv));
preKeys = cell(nPre, 1);
for i = 1:nPre
    preKeys{i} = keyOf(prePP{i}(:, 1));
end
lut = containers.Map();
for i = 1:nPre
    k = preKeys{i};
    if isKey(lut, k)
        lut(k) = [lut(k), i];
    else
        lut(k) = i;
    end
end

rows = cell(nPost, 1);
nUnique = 0; nAmb = 0; nUnmatched = 0; nGeomFail = 0; newOvrTotal = 0;
for j = 1:nPost
    P = postPP{j};
    idv = P(:, 1);
    k = keyOf(idv);
    cands = [];
    if isKey(lut, k); cands = lut(k); end
    hits = [];
    for c = cands
        if isequal(sort(prePP{c}(:, 1)), sort(idv))
            hits(end+1) = c; %#ok<AGROW>
        end
    end
    r = struct('postIdx', j, 'preIdx', NaN, 'matchType', {{'unmatched'}}, ...
        'nPP', size(P, 1), 'scale', NaN, 'maxDiamDev', NaN, 'maxCoordDev', NaN, ...
        'maxPairDistDev', NaN, 'newOverlapPairs', NaN, 'col6Equal', false, ...
        'rowOrderSame', false);
    if numel(hits) == 1
        r.preIdx = hits;
        r.matchType = {'id-exact'};
        nUnique = nUnique + 1;
        Q = prePP{hits};
        % Align rows by primary ID before comparing particle geometry.
        r.rowOrderSame = isequal(Q(:, 1), idv);
        if ~r.rowOrderSame
            [~, la] = ismember(idv, Q(:, 1));
            Q = Q(la, :);
        end
        ratio = P(:, 2) ./ Q(:, 2);
        sc = median(ratio);
        r.scale = sc;
        r.maxDiamDev = max(abs(ratio - sc)) / sc;
        w = Q(:, 2).^3;
        comQ = sum(w .* Q(:, 3:5), 1) ./ sum(w);
        w2 = P(:, 2).^3;
        comP = sum(w2 .* P(:, 3:5), 1) ./ sum(w2);
        cQ = (Q(:, 3:5) - comQ) .* sc;
        cP = P(:, 3:5) - comP;
        dref = sc * max(sqrt(sum((Q(:, 3:5) - comQ).^2, 2)) + Q(:, 2) ./ 2);
        if dref > 0
            r.maxCoordDev = max(sqrt(sum((cP - cQ).^2, 2))) / dref;
        else
            r.maxCoordDev = 0;
        end
        % Compare sampled pair distances after applying the diameter scale.
        n = size(P, 1);
        if n >= 2
            samp = unique(round(linspace(1, n, min(n, 120))));
            DP = pair_dists(P(samp, 3:5));
            DQ = pair_dists(Q(samp, 3:5)) .* sc;
            base = max(DQ, sc * min(Q(:, 2)));
            r.maxPairDistDev = max(abs(DP - DQ) ./ base);
        else
            r.maxPairDistDev = 0;
        end
        % Count post-scatter overlaps absent from the matched source aggregate.
        r.newOverlapPairs = count_new_overlap(Q, P);
        newOvrTotal = newOvrTotal + r.newOverlapPairs;
        r.col6Equal = isequal(Q(:, 6), P(:, 6));
        if max([r.maxDiamDev, r.maxCoordDev, r.maxPairDistDev]) > 1e-9
            r.matchType = {'id-exact-geomdev'};
            nGeomFail = nGeomFail + 1;
        end
    elseif numel(hits) > 1
        r.matchType = {'ambiguous'};
        r.preIdx = hits(1);
        nAmb = nAmb + 1;
    else
        nUnmatched = nUnmatched + 1;
    end
    rows{j} = struct2table(r);
end
T = vertcat(rows{:});
M = struct();
M.performed = true;
M.preFile = opts.PreScatterFile;
M.nPre = nPre;
M.nPost = nPost;
M.uniquelyMatched = nUnique;
M.ambiguous = nAmb;
M.unmatched = nUnmatched;
M.geometryDeviations = nGeomFail;
M.newOverlapPairsFromScatter = newOvrTotal;
M.maxDiamDev = max(T.maxDiamDev, [], 'omitnan');
M.maxCoordDev = max(T.maxCoordDev, [], 'omitnan');
M.maxPairDistDev = max(T.maxPairDistDev, [], 'omitnan');
M.allCol6Equal = all(T.col6Equal(strcmp(T.matchType, 'id-exact')));
M.allRowOrderSame = all(T.rowOrderSame(strcmp(T.matchType, 'id-exact')));
M.table = T;
logv(opts, ['  matched %d/%d uniquely (%d ambiguous, %d unmatched, ', ...
    '%d geometry deviations, %d new overlap pairs from scaling)\n'], ...
    nUnique, nPost, nAmb, nUnmatched, nGeomFail, newOvrTotal);
end

function D = pair_dists(xyz)
%PAIR_DISTS Return Euclidean distances for all unique row pairs.
n = size(xyz, 1);
[ii, jj] = find(triu(true(n), 1));
D = sqrt(sum((xyz(ii, :) - xyz(jj, :)).^2, 2));
end

function nNew = count_new_overlap(Q, P)
%COUNT_NEW_OVERLAP Count post-state overlap pairs absent from the pre-state.
%   Relative penetration greater than 1e-12 defines membership in each set.
sQ = overlap_set(Q);
sP = overlap_set(P);
nNew = size(setdiff(sP, sQ, 'rows'), 1);
end

function s = overlap_set(P)
%OVERLAP_SET Return row-index pairs that exceed the numerical threshold.
xyz = P(:, 3:5); d = P(:, 2);
cand = candidate_pairs_self(xyz, d);
res = exact_pairs(xyz, d, xyz, d, cand);
if isempty(res); s = zeros(0, 2); return; end
relp = -res(:, 5) ./ res(:, 4);
s = sortrows(res(relp > 1e-12, 1:2));
end

% ======================================================================= %
% Post-scatter -> parsdata(1) rigid-motion check
% ======================================================================= %

function rigid = post_to_initial_check(postPP, pd1PP, opts)
%POST_TO_INITIAL_CHECK Test whether LD2 initialization preserves geometry.
%   Corresponding cells must retain row count, IDs, diameters, and labels.
%   Coordinate differences must be constant within each aggregate, allowing
%   translation while rejecting rotation, scaling, or internal deformation.
logv(opts, '\n--- Post-scatter -> parsdata(1) rigid-transform check ---\n');
rigid = struct();
n = numel(postPP);
rigid.nAgg = n;
rigid.cellCountEqual = numel(pd1PP) == n;
if ~rigid.cellCountEqual
    rigid.note = 'cell counts differ; per-cell comparison skipped';
    return
end
badDiam = 0; badId = 0; badLbl = 0; badN = 0;
maxResid = 0; maxT = 0;
for i = 1:n
    A = postPP{i}; B = pd1PP{i};
    if size(A, 1) ~= size(B, 1); badN = badN + 1; continue; end
    if ~isequal(A(:, 1), B(:, 1)); badId = badId + 1; end
    if ~isequal(A(:, 2), B(:, 2)); badDiam = badDiam + 1; end
    if ~isequal(A(:, 6), B(:, 6)); badLbl = badLbl + 1; end
    dr = B(:, 3:5) - A(:, 3:5);
    t = mean(dr, 1);
    resid = max(sqrt(sum((dr - t).^2, 2)));
    maxResid = max(maxResid, resid);
    maxT = max(maxT, norm(t));
end
rigid.cellsWithSizeMismatch = badN;
rigid.cellsWithIdMismatch = badId;
rigid.cellsWithDiamMismatch = badDiam;
rigid.cellsWithLabelMismatch = badLbl;
rigid.maxTranslationResidual_m = maxResid;
rigid.maxTranslationNorm_m = maxT;
rigid.pureRigidTranslation = (badN + badId + badDiam + badLbl == 0) && ...
    maxResid < 1e-15; % Keep tolerance well below a typical primary radius.
logv(opts, ['  size/id/diam/label mismatches: %d/%d/%d/%d, ', ...
    'max translation residual %.3g m\n'], badN, badId, badDiam, badLbl, maxResid);
end

% ======================================================================= %
% Timeline and temporal tracing
% ======================================================================= %

function T = build_timeline(stages, stageStats, runMeta, opts) %#ok<INUSD>
%BUILD_TIMELINE Assemble stage order, population size, and overlap severity.
rows = cell(numel(stages), 1);
for s = 1:numel(stages)
    st = stageStats{s};
    r = struct();
    r.stage = {st.popId};
    r.kIter = st.kIter;
    r.nAgg = st.nAgg;
    r.pairsMinThr = st.pairsAtThr(1);
    r.pairsMaterial = st.pairsAtThr(min(3, numel(st.pairsAtThr)));
    r.maxRelPen = st.maxRelPen;
    rows{s} = struct2table(r);
end
T = vertcat(rows{:});
end

function trace = trace_pairs(stages, pairRows, idmaps, seedMap, maxId, opts)
%TRACE_PAIRS Track overlap-pair histories across ordered LD2 populations.
%   Histories are keyed by sorted primary-ID pairs and include the
%   post-scatter population. Duplicate primary IDs make the affected keys
%   ambiguous; exact geometric stage counts remain independent of this map.
ld2Mask = ~strcmp({stages.popId}, 'pre_scatter');
sIdx = find(ld2Mask);
keys = zeros(0, 2);
for s = sIdx
    P = pairRows{s};
    if isempty(P); continue; end
    k = sort([P.id1, P.id2], 2);
    keys = [keys; k]; %#ok<AGROW>
end
keys = unique(keys, 'rows');
nK = size(keys, 1);
nS = numel(sIdx);
relPen = nan(nK, nS);
coMember = false(nK, nS);
for c = 1:nS
    s = sIdx(c);
    P = pairRows{s};
    if ~isempty(P)
        k = sort([P.id1, P.id2], 2);
        [tf, loc] = ismember(k, keys, 'rows');
        relPen(loc(tf), c) = max(relPen(loc(tf), c), P.relPen(tf), 'omitnan');
        % Retain the greatest penetration when a key occurs more than once.
        for q = find(tf)'
            relPen(loc(q), c) = max([relPen(loc(q), c), P.relPen(q)], [], 'omitnan');
        end
    end
    map = idmaps{s};
    if ~isempty(map)
        [t1, l1] = ismember(keys(:, 1), map.id);
        [t2, l2] = ismember(keys(:, 2), map.id);
        ok = t1 & t2;
        coMember(ok, c) = map.agg(l1(ok)) == map.agg(l2(ok));
    end
end
firstCo = zeros(nK, 1);
firstOv = zeros(nK, 1);
for q = 1:nK
    f = find(coMember(q, :), 1);
    if ~isempty(f); firstCo(q) = f; end
    f = find(relPen(q, :) > 0, 1);
    if ~isempty(f); firstOv(q) = f; end
end
trace = struct();
trace.stageOrder = {stages(sIdx).popId};
trace.stageK = [stages(sIdx).kIter];
trace.pairKeys = keys;
trace.relPenByStage = relPen;
trace.coMemberByStage = coMember;
trace.firstCoMemberStage = firstCo;
trace.firstOverlapStage = firstOv;
% Measure the stage lag between first co-membership and first overlap.
born = zeros(nK, 1);
for q = 1:nK
    if firstOv(q) > 0 && firstCo(q) > 0
        born(q) = firstOv(q) - firstCo(q); % Zero denotes overlap at first co-membership.
    end
end
trace.overlapVsCoMemberLag = born;
% Attach source-seed identities when both primary IDs are unambiguous.
if ~isempty(seedMap)
    s1 = zeros(nK, 1); s2 = zeros(nK, 1);
    ok1 = keys(:, 1) >= 1 & keys(:, 1) <= maxId;
    ok2 = keys(:, 2) >= 1 & keys(:, 2) <= maxId;
    s1(ok1) = double(seedMap(keys(ok1, 1)));
    s2(ok2) = double(seedMap(keys(ok2, 2)));
    trace.seed1 = s1;
    trace.seed2 = s2;
    trace.crossSeedFrac = mean(s1 > 0 & s2 > 0 & s1 ~= s2);
end
logv(opts, '\nTemporal trace: %d distinct overlapping ID pairs across LD2 stages\n', nK);
end

% ======================================================================= %
% Deterministic property tests for repository collision functions
% ======================================================================= %

function A = test_connect_invariant(postPP, rs, opts)
%TEST_CONNECT_INVARIANT Evaluate the geometric postcondition of COL.CONNECT.
%   Stored aggregates are placed at prescribed positive and negative gaps
%   along deterministic directions. After capture, the test measures the
%   minimum residual gap and verifies that input arrays were not modified.
logv(opts, '\n--- Property test A: COL.CONNECT invariant ---\n');
A = struct('trials', [], 'summary', '');
if isempty(postPP)
    A.summary = 'skipped: post-scatter population unavailable';
    return
end
sizes = cellfun(@(x) size(x, 1), postPP);
targets = [8, 20, 50, 120, 300];
pickA = zeros(size(targets));
for q = 1:numel(targets)
    [~, pickA(q)] = min(abs(sizes - targets(q)));
end
pickB = zeros(size(targets));
for q = 1:numel(targets)
    [~, pickB(q)] = min(abs(sizes - targets(end + 1 - q)) + (1:numel(sizes))' * 1e-9);
end
dirs = [eye(3); [1 1 0] / sqrt(2); [1 1 1] / sqrt(3)];
for q = 1:4
    v = randn(rs, 1, 3);
    dirs = [dirs; v ./ norm(v)]; %#ok<AGROW>
end
rows = [];
for q = 1:numel(targets)
    P1 = center_pp(postPP{pickA(q)});
    P2 = center_pp(postPP{pickB(q)});
    dref = (max(P1(:, 2)) + max(P2(:, 2))) / 2;
    gaps = [0.8, 0.3, 0.05, 0, -0.25, -0.75, -1.5] .* dref;
    for u = 1:size(dirs, 1)
        for g = gaps
            [P2p, gAch] = place_at_gap(P1, P2, dirs(u, :), g);
            if isnan(gAch); continue; end
            in1 = P1; in2 = P2p;
            try
                [o1, o2, colstat] = COL.CONNECT(in1, in2);
            catch err
                r = base_rowA(q, u, g, gAch);
                r.error = {err.identifier};
                rows = [rows; struct2table(r)]; %#ok<AGROW>
                continue
            end
            r = base_rowA(q, u, g, gAch);
            r.error = {''};
            r.colstat = colstat;
            if colstat
                gAfter = min_cross_gap(o1, o2);
                r.minGapAfter = gAfter;
                r.selectedPairContact = abs(min_cross_gap_abs(o1, o2)) < 1e-9 * dref;
                r.residualOverlap = gAfter < -1e-12 * dref;
                r.relResidual = max(0, -gAfter) / dref;
            end
            r.inputsUntouched = isequal(in1, P1) && isequal(in2, P2p);
            rows = [rows; struct2table(r)]; %#ok<AGROW>
        end
    end
end
A.trials = rows;
cap = rows(rows.colstat == 1, :);
A.nCaptured = height(cap);
A.nResidualOverlap = nnz(cap.residualOverlap);
A.nResidualFromNonneg = nnz(cap.residualOverlap & cap.gTarget >= 0);
A.nResidualFromNeg = nnz(cap.residualOverlap & cap.gTarget < 0);
A.maxRelResidual = max([0; cap.relResidual]);
A.summary = sprintf(['%d captures: %d left residual overlap ', ...
    '(%d started non-overlapping, %d started overlapping); max residual %.3g x d_ref'], ...
    A.nCaptured, A.nResidualOverlap, A.nResidualFromNonneg, ...
    A.nResidualFromNeg, A.maxRelResidual);
logv(opts, '  %s\n', A.summary);
end

function r = base_rowA(q, u, g, gAch)
%BASE_ROWA Initialize one COL.CONNECT property-test result record.
r = struct('aggPair', q, 'dirIdx', u, 'gTarget', g, 'gAchieved', gAch, ...
    'colstat', 0, 'minGapAfter', NaN, 'selectedPairContact', false, ...
    'residualOverlap', false, 'relResidual', 0, 'inputsUntouched', false);
end

function P = center_pp(P)
%CENTER_PP Translate a primary-particle matrix to its mass-weighted center.
w = P(:, 2).^3;
com = sum(w .* P(:, 3:5), 1) ./ sum(w);
P(:, 3:5) = P(:, 3:5) - com;
end

function g = min_cross_gap(P1, P2)
%MIN_CROSS_GAP Return the minimum signed surface gap between two aggregates.
%   Negative values denote penetration and positive values denote separation.
g = inf;
blk = 256;
for i0 = 1:blk:size(P1, 1)
    i1 = min(i0 + blk - 1, size(P1, 1));
    ii = i0:i1;
    d = sqrt((P1(ii, 3) - P2(:, 3)').^2 + (P1(ii, 4) - P2(:, 4)').^2 + ...
        (P1(ii, 5) - P2(:, 5)').^2) - (P1(ii, 2) + P2(:, 2)') ./ 2;
    g = min(g, min(d(:)));
end
end

function g = min_cross_gap_abs(P1, P2)
%MIN_CROSS_GAP_ABS Return the smallest absolute cross-pair surface gap.
%   A value near zero indicates that at least one pair is at contact.
g = inf;
blk = 256;
for i0 = 1:blk:size(P1, 1)
    i1 = min(i0 + blk - 1, size(P1, 1));
    ii = i0:i1;
    d = sqrt((P1(ii, 3) - P2(:, 3)').^2 + (P1(ii, 4) - P2(:, 4)').^2 + ...
        (P1(ii, 5) - P2(:, 5)').^2) - (P1(ii, 2) + P2(:, 2)') ./ 2;
    g = min(g, min(abs(d(:))));
end
end

function [P2p, gAch] = place_at_gap(P1, P2, u, gTarget)
%PLACE_AT_GAP Translate centered P2 to a prescribed minimum cross-pair gap.
%   Bisection along direction u determines the translation. A negative target
%   deliberately creates an interpenetrating input for invariant testing.
R1 = max(sqrt(sum(P1(:, 3:5).^2, 2)) + P1(:, 2) ./ 2);
R2 = max(sqrt(sum(P2(:, 3:5).^2, 2)) + P2(:, 2) ./ 2);
tHi = R1 + R2 + abs(gTarget) + max(P1(:, 2));
tLo = 0;
gfun = @(t) min_cross_gap(P1, shift_pp(P2, t .* u));
gHi = gfun(tHi);
gLo = gfun(tLo);
if ~(gHi >= gTarget && gLo <= gTarget)
    P2p = []; gAch = NaN;
    return
end
for it = 1:80
    tm = (tHi + tLo) / 2;
    gm = gfun(tm);
    if gm > gTarget
        tHi = tm;
    else
        tLo = tm;
    end
end
P2p = shift_pp(P2, tHi .* u);
gAch = gfun(tHi);
end

function P = shift_pp(P, dr)
%SHIFT_PP Apply a Cartesian translation to primary-particle coordinates.
P(:, 3:5) = P(:, 3:5) + dr;
end

function B = test_duplicate_ids(postPP, opts)
%TEST_DUPLICATE_IDS Characterize COL.CONNECT behavior with repeated IDs.
%   Separate cases place a duplicate on the selected contact particle, on an
%   unrelated particle, and in the second aggregate.
logv(opts, '\n--- Property test B: duplicate primary IDs in COL.CONNECT ---\n');
B = struct();
% Construct a controlled fixture from a small stored aggregate.
sizes = cellfun(@(x) size(x, 1), postPP);
[~, k] = min(abs(sizes - 6));
P1 = center_pp(postPP{k});
P2 = center_pp(postPP{k});
P2(:, 1) = P2(:, 1) + 1e9; % Ensure that the two input ID sets are initially distinct.
% Place the clone inside the capture window with a positive surface gap.
[P2, gAch] = place_at_gap(P1, P2, [1, 0, 0], 0.3 * max(P1(:, 2)));
B.placementGap = gAch;
if isempty(P2)
    B.summary = 'skipped: placement failed';
    return
end
% Identify the rows that COL.CONNECT selects as the nearest cross pair.
[i1, i2] = nearest_cross_rows(P1, P2);
% Case B1: Duplicate the selected primary ID within the first aggregate.
P1dup = P1;
other = mod(i1, size(P1, 1)) + 1;
P1dup(other, 1) = P1dup(i1, 1);
B.case1_dupSelectedId = run_connect_capture(P1dup, P2);
% Case B2: Duplicate an unrelated primary ID within the first aggregate.
P1dup2 = P1;
oa = mod(i1, size(P1, 1)) + 1;
ob = mod(i1 + 1, size(P1, 1)) + 1;
if oa ~= i1 && ob ~= i1 && oa ~= ob
    P1dup2(ob, 1) = P1dup2(oa, 1);
end
B.case2_dupOtherId = run_connect_capture(P1dup2, P2);
% Case B3: Duplicate the selected primary ID within the second aggregate.
P2dup = P2;
o2 = mod(i2, size(P2, 1)) + 1;
P2dup(o2, 1) = P2dup(i2, 1);
B.case3_dupSelectedIdIn2 = run_connect_capture(P1, P2dup);
B.summary = sprintf('dup selected ID in pp1: %s | dup other ID: %s | dup selected ID in pp2: %s', ...
    B.case1_dupSelectedId.outcome{1}, B.case2_dupOtherId.outcome{1}, ...
    B.case3_dupSelectedIdIn2.outcome{1});
logv(opts, '  %s\n', B.summary);
end

function [i1, i2] = nearest_cross_rows(P1, P2)
%NEAREST_CROSS_ROWS Locate the cross pair with minimum signed surface gap.
best = inf; i1 = 1; i2 = 1;
for a = 1:size(P1, 1)
    d = sqrt(sum((P2(:, 3:5) - P1(a, 3:5)).^2, 2)) - (P1(a, 2) + P2(:, 2)) ./ 2;
    [m, b] = min(d);
    if m < best
        best = m; i1 = a; i2 = b;
    end
end
end

function out = run_connect_capture(P1, P2)
%RUN_CONNECT_CAPTURE Execute COL.CONNECT and normalize its test outcome.
out = struct('outcome', {{''}}, 'errId', {{''}}, 'colstat', NaN, 'minGapAfter', NaN);
try
    [o1, o2, colstat] = COL.CONNECT(P1, P2);
    out.colstat = colstat;
    if colstat
        out.minGapAfter = min_cross_gap(o1, o2);
        out.outcome = {'ok'};
    else
        out.outcome = {'no-capture'};
    end
catch err
    out.outcome = {'error'};
    out.errId = {err.identifier};
end
end

function C = test_chained_grow(opts)
%TEST_CHAINED_GROW Exercise sequential collision and index-remapping cases.
%   The fixtures isolate ordinary chaining, residual overlap after a merge,
%   candidate remapping, reversed-index mapping, and a complete GROW sequence
%   affected by an incorrect remap.
logv(opts, '\n--- Property test C: chained COL.GROW collisions ---\n');
C = struct();
og = struct('indupdate', 'off', 'col', 'agg');

% Scenario 1: Three-sphere chain with overlap in pairs A-B and B-C.
pp = {mk_pp(1, 1, [0, 0, 0], 11); mk_pp(2, 1, [0.8, 0, 0], 22); ...
    mk_pp(3, 1, [1.6, 0, 0], 33)};
C.scenario1 = run_grow_fixture(pp, og);

% Scenario 2: A merged pair {A,B} and particle C overlap at two locations.
% COL.CONNECT resolves only the cross pair with the greatest penetration.
ppAB = [mk_pp(1, 1, [0, 0, 0], 11); mk_pp(2, 1, [1, 0, 0], 11)];
ppC = mk_pp(3, 1, [0.55, 0.6, 0], 22);
C.scenario2 = run_grow_fixture({ppAB; ppC}, og);

% Scenario 3: Four-aggregate chain with candidates (1,2), (2,3), and (3,4).
pp = {mk_pp(1, 1, [0, 0, 0], 1); mk_pp(2, 1, [0.9, 0, 0], 2); ...
    mk_pp(3, 1, [1.8, 0, 0], 3); mk_pp(4, 1, [2.7, 0, 0], 4)};
C.scenario3 = run_grow_fixture(pp, og);

% Scenario 4: Validate COL.UNITE's returned map for a reversed pair i1 > i2.
% An earlier merge can remap a pending pair (a,b) to (a,b') with a > b'. The
% returned map must remain consistent with the physical cell deletion.
C.scenario4 = test_unite_reversed_pair();

% Scenario 5: Exercise the reversed-pair map through a complete GROW call.
% Candidates (1,5), (3,5), and (3,6) cause pending pair (3,5) to become the
% reversed pair (3,1) after UNITE(1,5). An incorrect i_list then remaps the
% final candidate (3,6) to different cells and leaves an overlap unresolved.
pp = {mk_pp(1, 1, [0, 0, 0], 1); mk_pp(2, 1, [10, 0, 0], 2); ...
    mk_pp(3, 1, [1.6, 0, 0], 3); mk_pp(4, 1, [20, 0, 0], 4); ...
    mk_pp(5, 1, [0.8, 0, 0], 5); mk_pp(6, 1, [2.4, 0, 0], 6)};
C.scenario5 = run_grow_fixture(pp, og);
C.scenario5.expectedAggIfMapCorrect = 3; % A6 should join the {A1,A5,A3} cluster.
C.scenario5.crossOverlapAfter = cross_overlap_report(C.scenario5.ppAfter);
C.scenario5.missedFlaggedMerge = C.scenario5.nAggAfter > 3 && ...
    C.scenario5.crossOverlapAfter.nPairs > 0;

C.summary = sprintf(['s1: %s | s2: %s (residual overlap %d, maxRelPen %.3g) | ', ...
    's3: %s | s4(reversed UNITE map): %s | s5(GROW misremap): nAgg=%d ', ...
    '(expected 3), unresolved flagged overlap pairs=%d'], ...
    C.scenario1.outcome{1}, C.scenario2.outcome{1}, ...
    C.scenario2.internalOverlapPairs, C.scenario2.maxRelPen, ...
    C.scenario3.outcome{1}, C.scenario4.summary, C.scenario5.nAggAfter, ...
    C.scenario5.crossOverlapAfter.nPairs);
logv(opts, '  %s\n', C.summary);
end

function S = test_unite_reversed_pair()
%TEST_UNITE_REVERSED_PAIR Validate COL.UNITE mapping for i_col = [3,1].
%   Unique primary IDs identify the physical destination of each original
%   cell, allowing direct comparison with the returned i_list mapping.
pp = cell(5, 1);
for i = 1:5
    pp{i} = mk_pp(100 + i, 1, [10 * i, 0, 0], i);
end
pars = struct('pp', {pp}, 'n', ones(5, 1), 'r', zeros(5, 3), ...
    'v', zeros(5, 3), 'm', ones(5, 1), 'dmax', ones(5, 1));
for i = 1:5
    pars.r(i, :) = pp{i}(1, 3:5);
end
S = struct('summary', '', 'iList', [], 'trueMap', [], 'corrupt', false, ...
    'wrongOldIndices', []);
try
    [pars2, iList] = COL.UNITE(pars, [3, 1], 'agg');
catch err
    S.summary = ['error: ', err.identifier];
    return
end
% Determine the surviving location of each original cell from its primary ID.
trueMap = nan(5, 1);
for i = 1:5
    for c = 1:numel(pars2.pp)
        if ismember(100 + i, pars2.pp{c}(:, 1))
            trueMap(i) = c;
            break
        end
    end
end
S.iList = iList;
S.trueMap = trueMap;
S.wrongOldIndices = find(iList(:, 2) ~= trueMap);
S.corrupt = ~isempty(S.wrongOldIndices);
if S.corrupt
    S.summary = sprintf('mapping mismatch (old indices [%s] misremapped)', ...
        num2str(S.wrongOldIndices'));
else
    S.summary = 'consistent';
end
end

function rep = cross_overlap_report(ppAfter)
%CROSS_OVERLAP_REPORT Measure remaining overlap between small test cells.
%   Property-test fixtures are intentionally small, so exhaustive aggregate
%   pair enumeration is appropriate for this check.
rep = struct('nPairs', 0, 'maxRelPen', 0);
nA = numel(ppAfter);
for a = 1:nA-1
    for b = a+1:nA
        P1 = ppAfter{a}; P2 = ppAfter{b};
        cutoff = (max(P1(:, 2)) + max(P2(:, 2))) / 2;
        cand = candidate_pairs_cross(P1(:, 3:5), P2(:, 3:5), cutoff);
        res = exact_pairs(P1(:, 3:5), P1(:, 2), P2(:, 3:5), P2(:, 2), cand);
        if isempty(res); continue; end
        relp = -res(:, 5) ./ res(:, 4);
        rep.nPairs = rep.nPairs + nnz(relp > 1e-9);
        rep.maxRelPen = max([rep.maxRelPen; relp]);
    end
end
end

function row = mk_pp(id, d, xyz, lbl)
%MK_PP Construct one synthetic primary-particle row.
row = [id, d, xyz, lbl];
end

function out = run_grow_fixture(pp, og)
%RUN_GROW_FIXTURE Build required aggregate metadata and execute COL.GROW.
%   The result records identity preservation, array alignment, and remaining
%   internal overlap after all reported collision processing.
pars = struct();
pars.pp = pp(:);
pars.n = cellfun(@(x) size(x, 1), pars.pp);
nA = numel(pars.pp);
pars.r = zeros(nA, 3);
pars.dmax = zeros(nA, 1);
for i = 1:nA
    w = pars.pp{i}(:, 2).^3;
    pars.r(i, :) = sum(w .* pars.pp{i}(:, 3:5), 1) ./ sum(w);
    pars.dmax(i) = max(2 .* sqrt(sum((pars.pp{i}(:, 3:5) - pars.r(i, :)).^2, 2)) ...
        + pars.pp{i}(:, 2));
end
pars.v = zeros(nA, 3);
pars.m = pars.n; % Use count-proportional placeholder masses for the fixture.
idsBefore = sort(cell2mat(cellfun(@(x) x(:, 1), pars.pp, 'UniformOutput', false)));
out = struct('outcome', {{'ok'}}, 'nAggAfter', NaN, 'idsPreserved', false, ...
    'aligned', false, 'internalOverlapPairs', 0, 'maxRelPen', 0, ...
    'errId', {{''}}, 'ppAfter', {{}});
try
    pars2 = COL.GROW(pars, og);
catch err
    out.outcome = {'error'};
    out.errId = {err.identifier};
    return
end
out.nAggAfter = numel(pars2.pp);
out.ppAfter = pars2.pp(:);
idsAfter = sort(cell2mat(cellfun(@(x) x(:, 1), pars2.pp(:), 'UniformOutput', false)));
out.idsPreserved = isequal(idsBefore, idsAfter);
out.aligned = numel(pars2.pp) == numel(pars2.n) && ...
    numel(pars2.pp) == size(pars2.r, 1) && numel(pars2.pp) == size(pars2.v, 1) && ...
    numel(pars2.pp) == numel(pars2.m) && ...
    all(cellfun(@(x) size(x, 1), pars2.pp(:)) == pars2.n(:));
for i = 1:numel(pars2.pp)
    P = pars2.pp{i};
    res = exact_pairs(P(:, 3:5), P(:, 2), P(:, 3:5), P(:, 2), ...
        candidate_pairs_self(P(:, 3:5), P(:, 2)));
    if isempty(res); continue; end
    relp = -res(:, 5) ./ res(:, 4);
    out.internalOverlapPairs = out.internalOverlapPairs + nnz(relp > 1e-9);
    out.maxRelPen = max([out.maxRelPen; relp]);
end
if ~out.idsPreserved || ~out.aligned
    out.outcome = {'index-corruption'};
end
end

function D = test_periodic(domSize, opts)
%TEST_PERIODIC Compare equivalent interior and cross-boundary collisions.
%   Both fixtures have the same minimum-image separation. Different merge
%   outcomes therefore identify inconsistent periodic collision handling.
logv(opts, '\n--- Property test D: periodic boundary handling in COL.GROW ---\n');
og = struct('indupdate', 'off', 'col', 'agg');
d = 1e-7;
L = domSize(:)';
% Place the interior pair at a center separation of 0.9d.
ppI = {mk_pp(1, d, [L(1)/2, L(2)/2, L(3)/2], 1); ...
    mk_pp(2, d, [L(1)/2 + 0.9*d, L(2)/2, L(3)/2], 2)};
D.interior = run_grow_fixture(ppI, og);
% Reproduce the same minimum-image separation across the positive x face.
ppX = {mk_pp(1, d, [0.45*d, L(2)/2, L(3)/2], 1); ...
    mk_pp(2, d, [L(1) - 0.45*d, L(2)/2, L(3)/2], 2)};
D.crossFace = run_grow_fixture(ppX, og);
D.interiorMerged = D.interior.nAggAfter == 1;
D.crossFaceMerged = D.crossFace.nAggAfter == 1;
D.missedPeriodicCollision = D.interiorMerged && ~D.crossFaceMerged;
D.summary = sprintf('interior merged: %d | cross-face merged: %d (minimum-image separation identical)', ...
    D.interiorMerged, D.crossFaceMerged);
logv(opts, '  %s\n', D.summary);
end

% ======================================================================= %
% Summaries, worst cases, rendering, saving
% ======================================================================= %

function T = summarize_stages(stageStats, opts)
%SUMMARIZE_STAGES Convert per-stage structures to a compact summary table.
%   Threshold-dependent columns report pair, aggregate, and primary-particle
%   counts. Remaining columns record severity and identifier classifications.
thr = opts.RelativeThresholds;
rows = cell(numel(stageStats), 1);
for s = 1:numel(stageStats)
    st = stageStats{s};
    r = struct();
    r.stage = {st.popId};
    r.kIter = st.kIter;
    r.file = {short_name(st.file)};
    r.nAgg = st.nAgg;
    r.nPP = st.nPP;
    r.minAggSize = st.minAggSize;
    r.maxAggSize = st.maxAggSize;
    for q = 1:numel(thr)
        r.(sprintf('pairs_rp%g', q)) = st.pairsAtThr(q);
        r.(sprintf('aggs_rp%g', q)) = st.aggsAtThr(q);
        r.(sprintf('pp_rp%g', q)) = st.ppAtThr(q);
    end
    r.fracAggAffected = st.aggsAtThr(1) / max(1, st.nAgg);
    r.maxAbsPen_m = st.maxAbsPen;
    r.maxRelPen = st.maxRelPen;
    r.medianRelPen = st.medianRelPen;
    r.p95RelPen = st.p95RelPen;
    r.p99RelPen = st.p99RelPen;
    r.sameLblPairs = st.sameLblPairs;
    r.crossLblPairs = st.crossLblPairs;
    r.dupIdPairs = st.dupIdPairs;
    r.minClearance_m = st.minClearance;
    r.contactPairs = st.contactPairs;
    r.roundoffOnlyPairs = st.roundoffOnlyPairs;
    rows{s} = struct2table(r);
end
T = vertcat(rows{:});
T.Properties.Description = sprintf('pairs_rp1..%d use thresholds [%s]', ...
    numel(thr), num2str(thr));
end

function s = short_name(f)
%SHORT_NAME Return a file name and extension without its parent path.
[~, n, e] = fileparts(f);
s = [n, e];
end

function T = cat_tables(c)
%CAT_TABLES Concatenate nonempty tables while preserving an empty result.
c = c(~cellfun(@(x) isempty(x) || (istable(x) && height(x) == 0), c));
if isempty(c)
    T = table();
else
    T = vertcat(c{:});
end
end

function W = collect_worst(worstPayload, opts)
%COLLECT_WORST Select distinct high-penetration pairs for rendering.
%   Sorted primary-ID pairs prevent the same physical pair from being
%   selected repeatedly when it persists across multiple saved stages.
all = [worstPayload{:}];
if isempty(all)
    W = struct([]);
    return
end
[~, ord] = sort([all.relPen], 'descend');
seen = containers.Map();
W = struct([]);
for q = ord
    P = all(q).pp;
    key = sprintf('%d_%d', sort([P(all(q).row1, 1), P(all(q).row2, 1)]));
    if isKey(seen, key); continue; end
    seen(key) = true;
    W = [W, all(q)]; %#ok<AGROW>
    if numel(W) >= opts.MaxRender; break; end
end
end

function figs = render_worst(W, opts)
%RENDER_WORST Create full-aggregate and local views of severe overlaps.
%   Subaggregate labels control color, while opacity and a connecting line
%   identify the selected primary-particle pair.
figs = gobjects(0);
[sx, sy, sz] = sphere(24);
for q = 1:min(numel(W), opts.MaxRender)
    w = W(q);
    P = w.pp;
    f = figure('Visible', 'off', 'Position', [50, 50, 1250, 620], 'Color', 'w');
    lbls = unique(P(:, 6));
    cmap = lines(max(numel(lbls), 2));
    for view_i = 1:2
        subplot(1, 2, view_i);
        hold on
        for i = 1:size(P, 1)
            ci = cmap(mod(find(lbls == P(i, 6), 1) - 1, size(cmap, 1)) + 1, :);
            focus = (i == w.row1) || (i == w.row2);
            if view_i == 2 && ~focus
                dpair = norm(P(w.row1, 3:5) - P(w.row2, 3:5));
                if norm(P(i, 3:5) - P(w.row1, 3:5)) > 4 * (dpair + P(w.row1, 2))
                    continue
                end
            end
            surf(P(i, 3) + sx .* P(i, 2) / 2, P(i, 4) + sy .* P(i, 2) / 2, ...
                P(i, 5) + sz .* P(i, 2) / 2, 'FaceColor', ci, ...
                'EdgeColor', 'none', 'FaceAlpha', tern(focus, 0.95, 0.12));
        end
        r1 = P(w.row1, 3:5); r2 = P(w.row2, 3:5);
        plot3([r1(1), r2(1)], [r1(2), r2(2)], [r1(3), r2(3)], 'k-', 'LineWidth', 2);
        % Mark the expected contact position of primary 2 along the center line.
        u = (r2 - r1) ./ max(norm(r2 - r1), eps);
        rc = r1 + u .* (P(w.row1, 2) + P(w.row2, 2)) ./ 2;
        plot3(rc(1), rc(2), rc(3), 'kx', 'MarkerSize', 14, 'LineWidth', 2);
        axis equal vis3d
        grid on
        view(3)
        camlight headlight
        lighting gouraud
        if view_i == 2
            mid = (r1 + r2) / 2;
            hw = 2.5 * (norm(r2 - r1) + max(P(w.row1, 2), P(w.row2, 2)));
            xlim(mid(1) + [-hw, hw]); ylim(mid(2) + [-hw, hw]); zlim(mid(3) + [-hw, hw]);
            title('Close-up');
        end
    end
    sgtitle(sprintf(['%s (k=%g) agg %d rows [%d,%d] ids [%d,%d] lbl [%d,%d] ', ...
        'clr=%.3g m relPen=%.3g'], w.stage, w.kIter, w.aggIdx, w.row1, w.row2, ...
        P(w.row1, 1), P(w.row2, 1), P(w.row1, 6), P(w.row2, 6), ...
        w.clearance, w.relPen), 'Interpreter', 'none');
    figs(end+1) = f; %#ok<AGROW>
end
end

function y = tern(c, a, b)
%TERN Select one of two scalar values from a logical condition.
if c; y = a; else; y = b; end
end

function outDir = save_outputs(report, figs, opts)
%SAVE_OUTPUTS Write tabular results, the report structure, and rendered cases.
%   Each invocation receives a timestamped directory so prior audits remain
%   unchanged.
stamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
outDir = fullfile(opts.OutputDir, ['audit__', stamp]);
if isfolder(outDir)
    outDir = [outDir, '_', num2str(randi(1e6))]; % Add a suffix on timestamp collision.
end
mkdir(outDir);
writetable(report.stageSummary, fullfile(outDir, 'stage_summary.csv'));
if istable(report.overlapPairs) && height(report.overlapPairs) > 0
    writetable(report.overlapPairs, fullfile(outDir, 'overlap_pairs.csv'));
end
if istable(report.indexIntegrity) && height(report.indexIntegrity) > 0
    writetable(report.indexIntegrity, fullfile(outDir, 'index_integrity.csv'));
end
if isfield(report.prePostScatterMatching, 'table')
    writetable(report.prePostScatterMatching.table, ...
        fullfile(outDir, 'pre_post_matching.csv'));
end
if istable(report.interaggregateSummary.perStage) && ...
        height(report.interaggregateSummary.perStage) > 0
    writetable(report.interaggregateSummary.perStage, ...
        fullfile(outDir, 'interaggregate_summary.csv'));
end
save(fullfile(outDir, 'overlap_audit_report.mat'), 'report', '-v7.3');
for q = 1:numel(figs)
    exportgraphics(figs(q), fullfile(outDir, sprintf('worst_case_%02d.png', q)), ...
        'Resolution', 140);
end
end

function txt = limitations_text(runMeta, report)
%LIMITATIONS_TEXT Document constraints on replay, timing, and provenance.
%   These statements accompany saved results so that audit conclusions are
%   interpreted within the available checkpoint and identifier history.
txt = { ...
    'MATLAB RNG state is not stored in the LD2 checkpoints (save() writes workspace variables only), so exact replay of any iteration interval is impossible; each resume restarted the default RNG stream.'
    sprintf('Only %d late checkpoints were retained (k >= 30001 of %d iterations); events between parsdata(1) (k=2) and the first checkpoint can only be bracketed, not pinpointed.', 4, runMeta.k_final)
    'parsdata(1) is overwritten during the first LD2 loop pass (the n_agg <= 1.0*n0 save condition is true at k=2), so the earliest saved LD2 state already includes one MARCH+PBC+GROW step; internal geometry is still rigid because MARCH/PBC only translate aggregates.'
    'The pre-scatter library aggregates have no common spatial frame, so interaggregate overlap is only evaluated for LD2-stage populations.'
    'Stage minimum-clearance values are computed over broad-phase candidate pairs (center distance below the per-aggregate maximum primary diameter); larger separations cannot be global minima when contacts exist.'
    'Pre->post provenance is established by exact column-1 ID-multiset matching plus geometric verification; if the merge library were regenerated with different ID offsets this matching would need revisiting.'
    'parsdata(2:5) save iterations are reconstructed from ensdata.n_agg threshold crossings (zero/unrecorded rows masked), not stored directly.'};
if isfield(runMeta, 'ensdataZeroRows') && ~isempty(runMeta.ensdataZeroRows)
    txt{end+1} = sprintf(['ensdata contains %d unrecorded iteration row(s) (k=[%s]) and cumulative-time ', ...
        'resets (k=[%s]): at least one additional mid-run interruption/resume occurred whose checkpoint ', ...
        'was not retained, so the ensdata.t axis is not a single continuous physical time.'], ...
        numel(runMeta.ensdataZeroRows), num2str(runMeta.ensdataZeroRows(:)'), ...
        num2str(runMeta.ensdataTimeResets(:)'));
end
if isfield(report, 'prePostScatterMatching') && ...
        isfield(report.prePostScatterMatching, 'unmatched') && ...
        report.prePostScatterMatching.unmatched > 0
    txt{end+1} = sprintf('%d post-scatter aggregates could not be matched to the candidate raw library; for these the file is only a candidate source.', ...
        report.prePostScatterMatching.unmatched);
end
end

function print_summary(report, opts)
%PRINT_SUMMARY Display concise stage counts and the earliest material overlap.
T = report.stageSummary;
thr = opts.RelativeThresholds;
iMat = find(thr <= 1e-3, 1, 'last');
if isempty(iMat); iMat = numel(thr); end
matCol = sprintf('pairs_rp%d', iMat);
aggCol = sprintf('aggs_rp%d', iMat);
fprintf('\n==================== AUDIT SUMMARY ====================\n');
fprintf('%-14s %6s %8s %10s %10s %11s\n', 'stage', 'nAgg', 'nPP', ...
    sprintf('pairs>=%g', thr(1)), sprintf('pairs>=%g', thr(iMat)), 'maxRelPen');
for s = 1:height(T)
    fprintf('%-14s %6d %8d %10d %10d %11.3g\n', T.stage{s}, T.nAgg(s), ...
        T.nPP(s), T.pairs_rp1(s), T.(matCol)(s), T.maxRelPen(s));
end
firstMat = find(T.(matCol) > 0, 1);
if isempty(firstMat)
    fprintf('\nNo material overlap (relPen >= %g) found at any saved stage.\n', thr(iMat));
else
    fprintf('\nEarliest saved stage with material overlap (relPen >= %g): %s (k=%g)\n', ...
        thr(iMat), T.stage{firstMat}, T.kIter(firstMat));
    fprintf('Affected aggregates there: %d of %d; worst relative penetration overall: %.3g\n', ...
        T.(aggCol)(firstMat), T.nAgg(firstMat), max(T.maxRelPen));
end
P = report.overlapPairs;
if istable(P) && height(P) > 0
    m = P.relPen >= thr(iMat);
    fprintf('Material pairs across all stages: %d (%d cross-label, %d same-label, %d involving duplicate IDs)\n', ...
        nnz(m), nnz(m & ~P.lblsEqual), nnz(m & P.lblsEqual), ...
        nnz(m & (P.id1DupInAgg | P.id2DupInAgg | P.idsEqual)));
    if any(strcmp('seed1', P.Properties.VariableNames))
        cs = m & P.seed1 > 0 & P.seed2 > 0 & P.seed1 ~= P.seed2;
        fprintf('Material pairs joining two different source seeds: %d\n', nnz(cs));
    end
end
fprintf('=======================================================\n');
end
