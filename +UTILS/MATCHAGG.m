function best_indices = MATCHAGG(parsdata, xfield, yfield, zfield, ...
    yweight, zweight, alpha, max_combinations, seed, x_max_raw)
% "MATCHAGG" selects one aggregate from each structure element in...
%   ...parsdata such that log(x) and y values are similar (weighted), and...
%   ...z values are as different as possible, penalizing raw x (npp) range.
% ----------------------------------------------------------------------- %
% Inputs:
% parsdata: N×1 struct array with fields for x, y, and z
% xfield: name of x field (e.g., 'npp'), log10 scaled
% yfield: name of y field (e.g., 'sigmapp'), linear
% zfield: name of z field (e.g., 'n_hyb'), integer
% yweight: weight for y similarity in [0, 1]
% zweight: weight for z diversity (nonnegative)
% alpha: penalty for npp (x) range across selected indices
% max_combinations: number of random index combinations to evaluate
% seed: random seed for reproducibility
% x_max_raw: maximum allowed raw value of xfield (e.g., npp)
% ----------------------------------------------------------------------- %
% Output:
% best_indices: 1×N vector of aggregate indices (one from each group)
% ----------------------------------------------------------------------- %

fprintf('>> MATCHAGG running with seed = %d and x_max_raw = %g\n', seed, x_max_raw)

N = numel(parsdata);

% gather all log10(x) and y values to normalize across all groups
all_logx = []; all_y = [];
for i = 1:N
    all_logx = [all_logx; log10(parsdata(i).(xfield)(:))];
    all_y = [all_y; parsdata(i).(yfield)(:)];
end
x_min = min(all_logx); x_max = max(all_logx);
y_min = min(all_y);    y_max = max(all_y);

% normalize and filter per group
for i = 1:N
    x_raw = parsdata(i).(xfield)(:);
    y_raw = parsdata(i).(yfield)(:);
    z_raw = parsdata(i).(zfield)(:);

    % apply x_max_raw filtering
    keep = x_raw <= x_max_raw;
    x_raw = x_raw(keep);
    y_raw = y_raw(keep);
    z_raw = z_raw(keep);

    if isempty(x_raw)
        error('MATCHAGG:GroupEmpty', ...
            'Group %d has no data after filtering with x_max_raw = %g', i, x_max_raw);
    end

    % normalize and store
    data(i).x_raw = x_raw;
    data(i).x_norm = (log10(x_raw) - x_min) / (x_max - x_min);
    data(i).y_raw = y_raw;
    data(i).y_norm = (y_raw - y_min) / (y_max - y_min);
    data(i).z = z_raw;
end

% initialize combinations
lens = arrayfun(@(d) length(d.x_raw), data);
rng(seed);
combos = zeros(max_combinations, N);
for k = 1:max_combinations
    for j = 1:N
        combos(k,j) = randi(lens(j));
    end
end

% loop through all combinations
best_score = -inf;
best_indices = nan(1, N);
for row = 1:max_combinations
    idx = combos(row, :);

    x_vals = arrayfun(@(i) data(i).x_norm(idx(i)), 1:N);
    y_vals = arrayfun(@(i) data(i).y_norm(idx(i)), 1:N);
    z_vals = arrayfun(@(i) data(i).z(idx(i)), 1:N);
    x_raw_vals = arrayfun(@(i) data(i).x_raw(idx(i)), 1:N);

    % compute similarity and diversity
    sim_cost = (1 - yweight) * std(x_vals) + yweight * std(y_vals);
    z_diff = sum(abs(z_vals' - z_vals), 'all') / 2;
    x_range = max(x_raw_vals) - min(x_raw_vals);

    % total score
    score = -sim_cost + zweight * z_diff - alpha * x_range;

    if score > best_score
        best_score = score;
        best_indices = idx;
    end
end

end
