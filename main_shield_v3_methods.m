%%% Post-process aggregates after second-stage Langevin dynamics to compute
%   the screening factor spp of each primary particle under three models:
%
%   METHODS (set `method` to 'opaque' | 'transparent' | 'semitransparent')
%
%   1) 'opaque'           (front-side, opaque):
%      A perimeter point on particle i is screened if at least one *other*
%      primary lies in front (z > z_i) and covers that point.
%
%   2) 'transparent'      (either-side, single layer sufficient):
%      Ignore z-ordering. A perimeter point is screened if it is covered by
%      at least one *other* primary in projection (total coverage ≥ 1).
%
%   3) 'semitransparent'  (projection, two-layer requirement):
%      Ignore z-ordering. A perimeter point is screened only if it is covered
%      by two or more *other* primaries in projection (total coverage ≥ 2).
%
%   spp in [0,1]: 0 = fully visible perimeter; 1 = fully screened.
%
%   Notes:
%   • This implementation assumes circular primary projections.
%   • Rotations below use intrinsic Euler angles ~ U[0,2π) for compatibility
%     with PAR.ROTATE. This is *not* uniform on SO(3). See the quaternion
%     stub below if you need unbiased random orientations later.

clc
clear
clf('reset')
close all

%% --- USER SETTINGS ---

% Input data
fdir_in  = 'D:\Hamed\CND\PhD\Publication\DLCA2\outputs\postLD2-01-Apr-2025_11-54-39_LD2-25NOV24';
fname_in = 'Post_LD2-25NOV24';

% Perimeter resolution
n_shield = 50;

% Screening method: 'opaque' | 'transparent' | 'semitransparent'
method = 'opaque';

% Options
opts.chunkSize = 512;   % tune to fit RAM/CPU; memory ~ chunkSize * n_shield booleans
opts.xyPrune   = true;  % fast reject if center distance > rk + rt (safe)

%% --- LOAD ---

S = load(fullfile(fdir_in, [fname_in '.mat']), 'parsdata');
parsdata = S.parsdata;

%% --- ROTATE RANDOMLY THEN COMPUTE SPP ---

n_shot = numel(parsdata);
theta  = linspace(0, 2*pi, n_shield);

for i = 1:n_shot
    n_agg = numel(parsdata(i).pp);
    parsdata(i).spp = cell(n_agg, 1);
    parsdata(i).spp_method = method;  % store selected method in this snapshot

    % Random aggregate rotations
    % NOTE: These are intrinsic Euler angles sampled uniformly in [0,2π),
    % which is *not* a uniform distribution on SO(3). This is kept to match
    % your existing PAR.ROTATE usage; consider unit quaternions later if you
    % need unbiased random orientations.
    %
    % --- Uniform random rotation (optional future improvement) ---
    % Shoemake's method for random unit quaternions (uniform on SO(3)):
    % u1 = rand(n_agg,1); u2 = rand(n_agg,1); u3 = rand(n_agg,1);
    % qx = sqrt(1-u1).*sin(2*pi*u2);
    % qy = sqrt(1-u1).*cos(2*pi*u2);
    % qz = sqrt(u1).*sin(2*pi*u3);
    % qw = sqrt(u1).*cos(2*pi*u3);
    % % Convert (qx,qy,qz,qw) to rotation matrices with quat2rotm() or an
    % % equivalent function, then rotate pp{j}(:,3:5) directly or adapt
    % % PAR.ROTATE to accept rotation matrices.
    %
    angs = 2*pi*rand(n_agg, 3);
    parsdata(i).pp = PAR.ROTATE(parsdata(i).pp, parsdata(i).npp, angs);

    fprintf('\nCalculating screening (%s) for snapshot %d of %d\n', method, i, n_shot);
    UTILS.TEXTBAR([0, n_agg]); UTILS.TEXTBAR([1, n_agg]);

    for j = 1:n_agg
        P   = parsdata(i).pp{j};      % columns: [.., d(2), x(3), y(4), z(5), ..]
        npp = parsdata(i).npp(j);

        % Preallocate spp for this aggregate
        spp = zeros(npp, 1);

        % Extract geometry
        d  = P(:,2);
        xc = P(:,3);
        yc = P(:,4);
        zc = P(:,5);
        r  = 0.5 * d;

        % Loop primaries
        for k = 1:npp
            rk = r(k);
            xk = xc(k) + rk * cos(theta);
            yk = yc(k) + rk * sin(theta);

            if npp <= 1
                spp(k) = 0;
                continue
            end

            % All other primaries (exclude self)
            tgt_all = setdiff(1:npp, k);

            % Optional XY pruning (safe, reduces work)
            if opts.xyPrune
                dx0  = xc(tgt_all) - xc(k);
                dy0  = yc(tgt_all) - yc(k);
                rtgt = r(tgt_all);
                keep = (dx0.^2 + dy0.^2) <= (rk + rtgt).^2;
                tgt_all = tgt_all(keep);
                if isempty(tgt_all)
                    spp(k) = 0;
                    continue
                end
            end

            % Counters per perimeter point (length n_shield)
            switch method
                case 'opaque'
                    n_plus  = zeros(1, n_shield, 'uint16'); % occluders with z > z_k

                case 'transparent'
                    n_all_1   = zeros(1, n_shield, 'uint16'); % total, ignore z

                case 'semitransparent'
                    n_all_2   = zeros(1, n_shield, 'uint16'); % total, ignore z

                otherwise
                    error('Unknown method: %s', method);
            end

            % Stream through targets in chunks (memory-safe)
            for c = 1:opts.chunkSize:numel(tgt_all)
                idx = tgt_all(c : min(c+opts.chunkSize-1, numel(tgt_all)));

                xt = xc(idx);
                yt = yc(idx);
                zt = zc(idx);
                rt = r(idx);

                % Vectorized inclusion test for this chunk: inside(j, m)
                dx = xk - xt;  % [numel(idx) x n_shield], implicit expansion
                dy = yk - yt;
                inside = (dx.^2 + dy.^2) <= (rt.^2);  % logical

                switch method
                    case 'opaque'
                        front = zt > zc(k);
                        if any(front)
                            n_plus = n_plus + uint16(sum(inside(front, :), 1));
                        end

                    case 'transparent'
                        % Ignore z; just total coverage by others
                        n_all_1 = n_all_1 + uint16(sum(inside, 1));

                    case 'semitransparent'
                        % Ignore z; coverage by more than one others
                        n_all_2 = n_all_2 + uint16(sum(inside, 1));
                end
            end

            % Final spp per method
            switch method
                case 'opaque'
                    % Screened if ≥1 front-side cover (viewer at +z)
                    spp(k) = mean(n_plus >= 1);

                case 'transparent'
                    % Screened if covered by ≥1 other (ignore z)
                    spp(k) = mean(n_all_1 >= 1);

                case 'semitransparent'
                    % Screened if covered by ≥2 others (ignore z)
                    spp(k) = mean(n_all_2 >= 2);
            end
        end

        parsdata(i).spp{j} = spp;
        UTILS.TEXTBAR([j, n_agg]);
    end
end

%% --- SAVE ---

ts      = char(datetime('now','Format','yyyy-MM-dd_HH-mm-ss'));
dir_out = fullfile('outputs', ['Shield_' ts '_' method]);
if ~isfolder(dir_out), mkdir(dir_out); end
save(fullfile(dir_out, ['Shield_' ts '_' method '.mat']), 'parsdata', 'n_shield', 'method', 'opts');
