clc
clear
clf('reset')
close all
warning('off')

%% initialize the script %%

% location of previously saved aggregate data
fdir_in = 'F:\DLCA2\outputs\postLD2-01-Apr-2025_11-54-39_LD2-25NOV24';
fname_in = 'Post_LD2-25NOV24';

varnames = {'parsdata'}; % avraiables to be imported

% load aggregate data
for i = 1 : numel(varnames)
    load(fullfile(fdir_in, strcat(fname_in, '.mat')), varnames{i});
end

ii0 = [1, 2, 3]; % snapshots of post-flame agglomeration to be plotted
n_agg_plt = length(ii0); % number of selected snapshots

% location of x-y-z axes in the aggregate plot
triad_pos = {[0.07, 0.71, 0.08, 0.08], [0.49, 0.71, 0.08, 0.08];
             [0.07, 0.42, 0.08, 0.08], [0.49, 0.42, 0.08, 0.08];
             [0.07, 0.11, 0.08, 0.08], [0.49, 0.11, 0.08, 0.08]};

% location of aggregate information on the plot
y_row_ttl = [0.96, 0.66, 0.37];
x_row_ttl = [0.2, 0.6];

% resolution for shielding calculation
n_shield = 50;

%% calculate shielding for all aggregates %%

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
              repmat(parsdata(i).pp{j}(:,2), 1, n_shield) .* ...
              cos(repmat(theta, parsdata(i).npp(j), 1));

        % y location
        y_c = parsdata(i).pp{j}(:,4) + ...
              repmat(parsdata(i).pp{j}(:,2), 1, n_shield) .* ...
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
                r_thresh2 = parsdata(i).pp{j}(tgt_idx,2).^2;
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

%% render selected aggregates - colorcode based on shielding %%

% identify indices of aggregates to be plotted
jj0 = UTILS.MATCHAGG(vertcat(parsdata(ii0)), 'npp', 'sigmapp', 'n_hyb', ...
    0.8, 80, 50.0, 10000, 38, 300);  % x_max_raw = 300

% display selected properties
for i = 1:n_agg_plt
    fprintf('Group %d: npp = %d, sigmapp = %.4f, n_hyb = %d\n', i, ...
        parsdata(ii0(i)).npp(jj0(i)), ...
        parsdata(ii0(i)).sigmapp(jj0(i)), ...
        parsdata(ii0(i)).n_hyb(jj0(i)));
end

% initialize figure
f2 = figure(2);
f2.Position = [100, 100, 600, 900];
set(f2, 'color', 'white')
tl2_tot = tiledlayout(n_agg_plt,2);
tl2_tot.TileSpacing = 'compact';
tl2_tot.Padding = 'compact';

% generate colormap
opts2.cc = 'on';
opts2.cm = colormap("abyss");
opts2.cm = flip(opts2.cm,1);

% placeholders for tiles
tl2 = cell(n_agg_plt, 2);
row_ttl = cell(3,1);

% render aggregates
for i = 1 : n_agg_plt

    % isometric view
    tl2{i,1} = nexttile(tl2_tot, 2*i-1);
    UTILS.PLOTPP_CONTINUOUS(parsdata(ii0(i)).pp{jj0(i)}(:,3),...
        parsdata(ii0(i)).pp{jj0(i)}(:,4),...
        parsdata(ii0(i)).pp{jj0(i)}(:,5),...
        parsdata(ii0(i)).pp{jj0(i)}(:,2),...
        parsdata(ii0(i)).spp{jj0(i)}, opts2);
    hold on
    UTILS.MAKEAX(f2, triad_pos{i,1}, '3d'); % x-y-z triad
    
    % x-y plane
    tl2{i,2} = nexttile(tl2_tot, 2*i);
    UTILS.PLOTPP_CONTINUOUS(parsdata(ii0(i)).pp{jj0(i)}(:,3),...
        parsdata(ii0(i)).pp{jj0(i)}(:,4),...
        parsdata(ii0(i)).pp{jj0(i)}(:,5),...
        parsdata(ii0(i)).pp{jj0(i)}(:,2),...
        parsdata(ii0(i)).spp{jj0(i)}, opts2);
    hold on
    view(2)
    UTILS.MAKEAX(f2, triad_pos{i,2}, 'xy');
    
    % print aggregate structural information
    row_ttl{i} = strcat('$n_\mathrm{pp}$ =', {' '},...
        num2str(parsdata(ii0(i)).npp(jj0(i))),...
        ', $\sigma_\mathrm{pp}$ =', {' '},...
        num2str(parsdata(ii0(i)).sigmapp(jj0(i)), '%.2f'),...
        ', $n_\mathrm{hyb}$ =', {' '},...
        num2str(parsdata(ii0(i)).n_hyb(jj0(i))));
    annotation('textbox', [x_row_ttl(1), y_row_ttl(i), x_row_ttl(2), 0.03],...
        'String', row_ttl{i}, 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'bottom', 'FontSize', 14, 'EdgeColor',...
        'none', 'Interpreter', 'latex');
    
end

% generate colorbar showing shielding values
cb2 = colorbar(tl2{1,1}, 'eastoutside');
cb2.Layout.Tile = 'south';
cb2.Label.String = 'Shielding ratio [-]';
cb2.FontSize = 12;
cb2.Label.FontSize = 18;
cb2.TickLabelInterpreter = 'latex';
cb2.Label.Interpreter = 'latex';
cb2.LineWidth = 1;
