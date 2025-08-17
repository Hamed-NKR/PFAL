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
n_shield = 64;

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

        parsdata(i).spp{j} = PAR.SHIELD(parsdata(i).pp{j}, n_shield,...
            method, opts); % screening factor for each individual aggregate

        UTILS.TEXTBAR([j, n_agg]); % update progress bar 

    end
    
end

%% --- SAVE ---

ts      = char(datetime('now','Format','yyyy-MM-dd_HH-mm-ss'));
dir_out = fullfile('outputs', ['Shield_' ts '_' method]);
if ~isfolder(dir_out), mkdir(dir_out); end
save(fullfile(dir_out, ['Shield_' ts '_' method '.mat']), 'parsdata', 'n_shield', 'method', 'opts');
