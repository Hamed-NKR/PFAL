%%% A script that post-processes aggregates after second-stage Langevin
%   dynamics to calculate the screening factor of each individual
%   primary particle within them. This is achieved through:
%   (a) randomly rotating aggregates, (b) discretizing the perimeters
%   of their primary-particle projections along the z axis, and
%   (c) checking whether the discretized perimeter points fall within
%   the projected area of at least two *other* primary particles when
%   observed from the z orientation.
%
%   Finally, if ≥50% of the discretized perimeter points on primary
%   particle [i] lie within the z-projection area of at least two other
%   primary particles in the same aggregate, primary particle [i] is
%   considered screened (i.e., “invisible”) for that random view.
%   The choice of z-axis view should not bias the calculation because
%   aggregates are randomly rotated first. A new field in the particle
%   data structure is generated to store the screening factor, denoted
%   by spp.
%
%   "pars.spp = 0" means the entire perimeter of the primary particle is
%   visible from a random projection, whereas "pars.spp = 1" means the
%   primary-particle perimeter is fully covered by other primary particles
%   in that view. Why require overlap by “two” other projections?
%   The goal is to mimic manual TEM sizing: because real primaries
%   often appear more elliptical than circular, once >50% of a primary’s
%   perimeter is overlapped by two or more neighbors, the perimeter is
%   visually ambiguous and typically skipped in manual sizing.
%
%   Notes:
%   • TEM images are projections: there is no nearer/farther occlusion
%     order along z. We therefore do not depth-filter; any overlapping
%     projection contributes to screening.
%   • This implementation assumes circular primary projections (radius = d/2).
%     If you later adopt ellipses, replace the circle inclusion test with
%     an ellipse implicit form after projection.

clc
clear
clf('reset')
close all
% warning('off')  % Keep warnings ON while developing; disable only if needed.

%% Initialize the script %%

% Location of previously saved aggregate data
fdir_in  = 'D:\Hamed\CND\PhD\Publication\DLCA2\outputs\postLD2-01-Apr-2025_11-54-39_LD2-25NOV24';
fname_in = 'Post_LD2-25NOV24';

% Load aggregate data
S = load(fullfile(fdir_in, [fname_in '.mat']), 'parsdata');
parsdata = S.parsdata;

% Resolution for perimeter discretization in screening calculations
n_shield = 50;

%% Calculate screening factor for all aggregates across post-flame snapshots %%

n_shot = numel(parsdata);          % total number of post-flame snapshots
n_agg  = zeros(n_shot,1);          % number of aggregates per snapshot

% Discretize perimeter angle for screening calculation
theta = linspace(0, 2*pi, n_shield);

for i = 1:n_shot

    % Number of aggregates in this post-flame snapshot
    n_agg(i) = numel(parsdata(i).pp);

    % Allocate a field for screening factors in each aggregate
    parsdata(i).spp = cell(n_agg(i),1);

    % Random aggregate rotations
    % NOTE: These are intrinsic Euler angles sampled uniformly in [0,2π),
    % which is *not* a uniform distribution on SO(3). This is kept to match
    % your existing PAR.ROTATE usage; consider unit quaternions later if you
    % need unbiased random orientations.
    %
    % --- Uniform random rotation (optional future improvement) ---
    % Shoemake's method for random unit quaternions (uniform on SO(3)):
    % u1 = rand(n_agg(i),1);
    % u2 = rand(n_agg(i),1);
    % u3 = rand(n_agg(i),1);
    % qx = sqrt(1-u1).*sin(2*pi*u2);
    % qy = sqrt(1-u1).*cos(2*pi*u2);
    % qz = sqrt(u1).*sin(2*pi*u3);
    % qw = sqrt(u1).*cos(2*pi*u3);
    % % Convert (qx,qy,qz,qw) to rotation matrices with quat2rotm() or an
    % % equivalent function, then rotate coordinates directly or adapt
    % % PAR.ROTATE to accept rotation matrices.
    %
    angs = 2*pi*rand(n_agg(i), 3);

    % Rotate aggregates using the angles above
    parsdata(i).pp = PAR.ROTATE(parsdata(i).pp, parsdata(i).npp, angs);

    disp(' ')
    fprintf('Calculating screening for post-flame snapshot %d of %d\n', i, n_shot);

    UTILS.TEXTBAR([0, n_agg(i)]);
    UTILS.TEXTBAR([1, n_agg(i)]); % start progress textbar

    for j = 1:n_agg(i)

        % Initialize screening factor for each primary particle in this aggregate
        npp = parsdata(i).npp(j);
        P   = parsdata(i).pp{j};  % columns assumed: [.., d(2), x(3), y(4), z(5), ..]
        parsdata(i).spp{j} = zeros(npp, 1);

        % Extract centers and radii once
        d  = P(:,2);
        xc = P(:,3);
        yc = P(:,4);
        r  = 0.5*d;

        % Choose a chunk size that fits comfortably in RAM
        % Memory per chunk ~ chunkSize * n_shield booleans + overhead
        chunkSize = 512;

        % Loop over primaries in this aggregate
        for k = 1:npp

            % Perimeter samples for particle k (compute once per k)
            rk = r(k);
            xk = xc(k) + rk * cos(theta);
            yk = yc(k) + rk * sin(theta);

            % Candidates = all other particles (no z filtering in TEM projection)
            if npp > 1
                all_idx = 1:npp;
                tgt_all = all_idx(all_idx ~= k).';
            else
                parsdata(i).spp{j}(k) = 0;
                continue
            end

            % Optional XY pruning to safely cut unnecessary tests:
            % If center distance > rk + rt, that target cannot cover any perimeter point.
            dx0  = xc(tgt_all) - xc(k);
            dy0  = yc(tgt_all) - yc(k);
            rtgt = r(tgt_all);
            keep = (dx0.^2 + dy0.^2) <= (rk + rtgt).^2;
            tgt_all = tgt_all(keep);

            % Running count of how many *other* primaries cover each perimeter point
            % Use uint16 for compactness (n_cover per point rarely exceeds a few dozen)
            n0_points = zeros(1, n_shield, 'uint16');

            % Stream through targets in chunks to keep memory usage bounded
            for c = 1:chunkSize:numel(tgt_all)
                idx = tgt_all(c : min(c+chunkSize-1, numel(tgt_all)));

                xt = xc(idx);
                yt = yc(idx);
                rt = r(idx);

                % Vectorized inclusion test for this chunk:
                % inside is [numel(idx) x n_shield] logical
                dx = xk - xt;      % implicit expansion
                dy = yk - yt;
                inside = (dx.^2 + dy.^2) <= (rt.^2);

                % Accumulate occluder counts per perimeter point (column-wise sum)
                n0_points = n0_points + uint16(sum(inside, 1));
            end

            % Screening fraction: fraction of perimeter points covered by ≥ 2 others
            parsdata(i).spp{j}(k) = mean(n0_points >= 2);
        end

        UTILS.TEXTBAR([j, n_agg(i)]); % update progress textbar
    end
end

%% Save updated aggregate data including screening factors %%

% Make an output directory and save only what’s needed
ts      = char(datetime('now','Format','yyyy-MM-dd_HH-mm-ss'));
dir_out = fullfile('outputs', ['Shield_' ts]);
if ~isfolder(dir_out)
    mkdir(dir_out);
end

save(fullfile(dir_out, ['Shield_' ts '.mat']), 'parsdata', 'n_shield');
