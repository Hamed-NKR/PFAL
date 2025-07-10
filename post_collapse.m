clc
clear
clf('reset')
close all
warning('off')

%% initialize %%

% load aggregate data
fdir_in = 'C:\Users\hmdni\Downloads\Colaps';
fname_in = 'Colaps';
varnames = {'pars_clps','pars_out'}; % varaiables to be imported
for i = 1 : numel(varnames)
    load(fullfile(fdir_in, strcat(fname_in, '.mat')), varnames{i});
end

pars_clps = {pars_out, pars_clps};

% assign resolution for shielding calculations
n_shield = 20;

% discretize orientation around primary particles for shielding calculation
theta = linspace(0, 2*pi, n_shield);

% total number of aggregates
% n_agg = length(pars_clps);
n_agg = 748;

%% calculate shielding factor %%

n_clps_shot = length(pars_clps);

for i = 1 : n_clps_shot 

    % allocate a field for shielding for each aggregate
    pars_clps{i}.spp = cell(n_agg,1);

    angs = 2 * pi * rand(n_agg, 3); % a uniform random set of three...
    % ...intrinsic Euler angles for aggregte rotation specific to...
    % ...individual aggregates

    % % rotate aggregates in a random direction given above
    % pars_clps{i}.pp = PAR.ROTATE(pars_clps{i}.pp, pars_clps{i}.npp, angs);

    disp(' ')
    fprintf('Post-flame snapshot %d of %d\n', i, n_clps_shot);

    UTILS.TEXTBAR([0, n_agg]);
    UTILS.TEXTBAR([1, n_agg]); % start progress textbar

    for j = 1 : n_agg

        % initialize shielding factor for individual primary particles
        pars_clps{i}.spp{j} = zeros(pars_clps{i}.n(j), 1);

        % calculate shielding factor (with spp: [0 1] where spp = 1...
        % ...meaning fully screened) for each primary particle in...
        % ...each aggregate

        % make primary particle pair indices
        kk = table2array(combinations(1 : pars_clps{i}.n(j), ...
            1 : pars_clps{i}.n(j)));

        % whether a primary particle is further from the viewer than...
        % ...another primary particle
        dz = pars_clps{i}.pp{j}(kk(:,1),5) < pars_clps{i}.pp{j}(kk(:,2),5);

        % locations of circumferential points on perimeter of each...
        % ...primary particle in x-y plane

        % x location of circumferential points
        x_c = pars_clps{i}.pp{j}(:,3) + ...
            repmat(pars_clps{i}.pp{j}(:,2)/2, 1, n_shield) .* ...
            cos(repmat(theta, pars_clps{i}.n(j), 1));

        % y location
        y_c = pars_clps{i}.pp{j}(:,4) + ...
            repmat(pars_clps{i}.pp{j}(:,2)/2, 1, n_shield) .* ...
            sin(repmat(theta, pars_clps{i}.n(j), 1));

        % set the maximum amount of shielding for each primary particle...
        % ...from all other primary particles as the final value...
        % ...assigned for shielding

        chunkSize = 5000; % adjust this based on available memory

        for k = 1 : pars_clps{i}.n(j)

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
                dx = x_c(src_idx,:) - pars_clps{i}.pp{j}(tgt_idx,3);
                dy = y_c(src_idx,:) - pars_clps{i}.pp{j}(tgt_idx,4);
                r2 = dx.^2 + dy.^2;

                % compare to radius of other particle (squared threshold)
                r_thresh2 = pars_clps{i}.pp{j}(tgt_idx,2).^2 / 4;
                r_thresh2 = repmat(r_thresh2, 1, size(dx,2));

                % whether circumferential points fall within projections...
                % ...of other primary particles
                dr_chunk = r2 < r_thresh2;

                % count the number of shielded circumferential points
                spp0(c:c_end) = sum(dr_chunk, 2);
            end

            % select maximum value as final shielding
            if ~isempty(spp0)
                pars_clps{i}.spp{j}(k) = max(spp0) / n_shield;
            end

        end

        UTILS.TEXTBAR([j, n_agg]); % update progress textbar

    end

end

%% post-process

% initialize figure
f1 = figure(1);
f1.Position = [50, 50, 1050, 550];
set(f1, 'color', 'white');

% initialize layout
tl1 = tiledlayout(1, 2);
tl1.TileSpacing = 'compact';
tl1.Padding = 'loose';

% assign colors for pots-flame snapshots
clr1 = colormap(hot);
cind1 = round(1 + (length(clr1) - 1) .* (0.05 : 0.7 / (5 - 1) : 0.75)');
clr1 = clr1(cind1,:);
clr1(end,:) = [236,230,61] / 255;

nexttile(1)

n_clps_shot = 2;

spp = cell(n_clps_shot, 1); % allocate variable for ensemble shielding factors
mu_spp = zeros(n_clps_shot, 1); % allocate means of shielding factor

% allocate probability density function (PDF) and evaluation points for...
    % ...ensemble shielding
f_spp = cell(n_clps_shot, 1);
xi_spp = cell(n_clps_shot, 1);

scale_spp = -0.15; % assign a scale for PFD to adjust curve extents

% allocate labels for x axis
xlbl1 = cell(n_clps_shot, 1);

npp_tot = zeros(n_clps_shot, 1); % allocate variable for the ensemble number...
    % ...of primary particles

hold on
for i = 1 : n_clps_shot

    % concatinate shielding factors across aggregates
    spp{i} = cat(1, pars_clps{i}.spp{:});    
    
    % compute mean of shielding factor for each snapshot
    mu_spp(i) = mean(spp{i});

    % % make a label for the distribution based on post-flame lifetime
    % if i == 1
    %     xlbl1{i} = num2str(pars_clps{i}.r_n_agg(1), '%.0f');
    % elseif ismember(i, [2,3])
    %     xlbl1{i} = num2str(pars_clps{i}.r_n_agg(1), '%.1f');
    % else
    %     xlbl1{i} = num2str(pars_clps{i}.r_n_agg(1), '%.2f');
    % end
    
    scatter(i-0.5, 1.2) % this is only to assign xticklabels
    
    if i == 1
        % print mean of each distribution
        text((i - 0.82), 1.03, ...
            sprintf('$\\langle{S}_\\mathrm{pp}\\rangle$ = %.2f', mu_spp(i)),...
            'Interpreter', 'latex', 'HorizontalAlignment', 'center',...       
            'FontSize', 16)
    else
        % print mean of each distribution
        text((i - 0.48), 1.03, num2str(mu_spp(i), '%.2f'),...
            'Interpreter', 'latex', 'HorizontalAlignment', 'center',...       
            'FontSize', 16)    
    end
    
    % total number of primary particles across all aggregates
    npp_tot(i) = length(spp{i});
    
    % evaluate a probability density function od shielding factor for...
        % ...the entire population of primary particles
    [f_spp{i}, xi_spp{i}] = ksdensity(spp{i});
    
    % generate proper data format for shading beneath the distribution
    y_fill = [scale_spp * f_spp{i} + (i - 0.5),...
        (i - 0.5) * ones(size(f_spp{i}))];
    x_fill = [xi_spp{i}, fliplr(xi_spp{i})];
    
    % % plot shielding distribution for initial snapshot (for comparison)
    % if i == 1
    %     % save initial distribution
    %     y_fill_0 = y_fill;
    %     x_fill_0 = x_fill;
    % else
    %     plot(scale_spp * f_spp{1} + (i - 0.5), xi_spp{1}, 'Color',...
    %         [clr1(1,:), 0.5], 'LineWidth', 0.5);
    %     fill(y_fill_0 + (i - 1), x_fill_0, clr1(1,:), ...
    %         'FaceAlpha', 0.1, ...
    %         'EdgeColor', 'none')
    % end

    % plot shielding distribution for current snapshot
    plot(scale_spp * f_spp{i} + (i - 0.5), xi_spp{i}, 'Color', clr1(i,:),...
        'LineWidth', 2);
    fill(y_fill, x_fill, clr1(i,:), ...
        'FaceAlpha', 0.3, ...
        'EdgeColor', 'none')

end

xticks((1 : n_clps_shot) - 0.5)  % specify tick positions for horizontal axis
% xticklabels(xlbl1)  % assign labels to ticks

% set plot appearances
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16,...
    'TickLength', [0.02 0.02])
xlim([-0.5 2])
ylim([0 1])
xlabel('$n_\mathrm{agg}/(n_\mathrm{agg})_2$ [-]', 'interpreter', 'latex',...
    'FontSize', 22)
ylabel('$S_\mathrm{pp}^\mathrm{(i)}$ [-]', 'interpreter', 'latex',...
    'FontSize', 22)

%% save updated aggregate data inlcuding shielding factor %%

% make a directory to save outputs
dir0_out = datestr(datetime('now'));
dir0_out = regexprep(dir0_out, ':', '-');
dir0_out = regexprep(dir0_out, ' ', '_');
dir_out = strcat('outputs\', 'Collapse_plus_Shield_', dir0_out, '\');
if ~isfolder(dir_out)
    mkdir(dir_out); % if it doesn't exist, create the directory
end

% save worksapce
save(strcat(dir_out, 'Collapse_plus_Shield_', dir0_out, '.mat'))
