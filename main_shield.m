%%% a script that post-process aggregates after second-stage Langevin...
    % ...dynamics simulation to calculate the shielding factor of their...
    % ...indidvidual primary particles. This is performed via...
    % ...(i) randomly rotating aggregates, (ii) descritizing the...
    % ...perimeters of their primary particle projections,...
    % ...(iii) deciding wether the descretized points fall below another...
    % ...primary particles when observed from a certain orientation. The...
    % ...observation view in this case is positive z axis. This specific...
    % ...view shouldn't make any bias becasue aggregates are first...
    % ...randomly rotated. A new field in particle data structure is...
    % ...generated to store shielding factor, named spp. spp=0 means...
    % ...all perimeter of primary particle is visible from the viewed...
    % ...projection, and spp=1 means primary particle is fully covered...
    % ...by other primary particles in that certain view.

clc
clear
clf('reset')
close all
warning('off')

%% initialize the script %%

% location of previously saved aggregate data
fdir_in = 'F:\DLCA2\outputs\postLD2-01-Apr-2025_11-54-39_LD2-25NOV24';
fname_in = 'Post_LD2-25NOV24';

varnames = {'parsdata'}; % varaiables to be imported

% load aggregate data
for i = 1 : numel(varnames)
    load(fullfile(fdir_in, strcat(fname_in, '.mat')), varnames{i});
end

% assign resolution for shielding calculations
n_shield = 50;

%% calculate shielding factor for all aggregates from various post-flame snapshots %%

n_shot = length(parsdata); % total number of post-flame snapshots
n_agg = zeros(n_shot,1); % allocate space to number of aggregate

% discretize orientation around primary particles for shielding calculation
theta = linspace(0, 2*pi, n_shield);

for i = 1 : n_shot

    % number of aggregates in each post-flame snapshot
    n_agg(i) = length(parsdata(i).pp);

    % allocate a field for shielding for each aggregate
    parsdata(i).spp = cell(n_agg(i),1);

    angs = 2 * pi * rand(n_agg(i), 3); % a uniform random set of three...
        % ...intrinsic Euler angles for aggregte rotation specific to...
        % ...individual aggregates
    
    % rotate aggregates in a random direction given above
    parsdata(i).pp = PAR.ROTATE(parsdata(i).pp, parsdata(i).npp, angs);

    disp(' ')
    fprintf('Post-flame snapshot %d of %d\n', i, n_shot);

    UTILS.TEXTBAR([0, n_agg(i)]);
    UTILS.TEXTBAR([1, n_agg(i)]); % start progress textbar

    for j = 1 : n_agg(i)

        % initialize shielding factor for individual primary particles
        parsdata(i).spp{j} = zeros(parsdata(i).npp(j), 1);

        % calculate shielding factor (with spp: [0 1] where spp = 1...
            % ...meaning fully screened) for each primary particle in...
            % ...each aggregate

        % make primary particle pair indices
        kk = table2array(combinations(1 : parsdata(i).npp(j), ...
                                      1 : parsdata(i).npp(j)));

        % whether a primary particle is further from the viewer than...
            % ...another primary particle
        dz = parsdata(i).pp{j}(kk(:,1),5) < parsdata(i).pp{j}(kk(:,2),5);

        % locations of circumferential points on perimeter of each...
            % ...primary particle in x-y plane

        % x location of circumferential points
        x_c = parsdata(i).pp{j}(:,3) + ...
              repmat(parsdata(i).pp{j}(:,2)/2, 1, n_shield) .* ...
              cos(repmat(theta, parsdata(i).npp(j), 1));

        % y location
        y_c = parsdata(i).pp{j}(:,4) + ...
              repmat(parsdata(i).pp{j}(:,2)/2, 1, n_shield) .* ...
              sin(repmat(theta, parsdata(i).npp(j), 1));

        % set the maximum amount of shielding for each primary particle...
            % ...from all other primary particles as the final value...
            % ...assigned for shielding

        chunkSize = 5000; % adjust this based on available memory

        for k = 1 : parsdata(i).npp(j)

            % find indices of pair-wise combination specific to each...
                % ...primary particle
            kkk = find((kk(:,1) == k) & dz); % choose only the ones...
                % ...with lower z

            % delete raw values of number of primary particles shielded...
                % ...from previous iteration 
            spp0 = zeros(1, length(kkk));

            % process in chunks to reduce memory usage
            for c = 1 : chunkSize : length(kkk)
                c_end = min(c + chunkSize - 1, length(kkk));
                idx_chunk = kkk(c:c_end);
                src_idx = kk(idx_chunk,1); % always equal to k
                tgt_idx = kk(idx_chunk,2);

                % circumferential point distances to other primary...
                    % ...particle centers
                dx = x_c(src_idx,:) - parsdata(i).pp{j}(tgt_idx,3);
                dy = y_c(src_idx,:) - parsdata(i).pp{j}(tgt_idx,4);
                r2 = dx.^2 + dy.^2;

                % compare to radius of other particle (squared threshold)
                r_thresh2 = parsdata(i).pp{j}(tgt_idx,2).^2 / 4;
                r_thresh2 = repmat(r_thresh2, 1, size(dx,2));

                % whether circumferential points fall within projections...
                    % ...of other primary particles
                dr_chunk = r2 < r_thresh2;

                % count the number of shielded circumferential points
                spp0(c:c_end) = sum(dr_chunk, 2);
            end

            % select maximum value as final shielding
            if ~isempty(spp0)
                parsdata(i).spp{j}(k) = max(spp0) / n_shield;
            end

        end

        UTILS.TEXTBAR([j, n_agg(i)]); % update progress textbar

    end

end

%% save updated aggregate data inlcuding shielding factor %%

% make a directory to save outputs
dir0_out = datestr(datetime('now'));
dir0_out = regexprep(dir0_out, ':', '-');
dir0_out = regexprep(dir0_out, ' ', '_');
dir_out = strcat('outputs\', 'Shield_', dir0_out, '\');
if ~isfolder(dir_out)
    mkdir(dir_out); % if it doesn't exist, create the directory
end

% save worksapce
save(strcat(dir_out, 'Shield_', dir0_out, '.mat'))
