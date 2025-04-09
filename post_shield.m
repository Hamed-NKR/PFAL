clc
clear
clf('reset')
close all
warning('off')

%% initialize the script %%

% location of previously saved aggregate data (this should include...
    % ...spp field in pars strcuture for primary particle shielding)
fdir_in = 'C:\Users\hmdnkr\Documents\GitHub\MCEM\outputs\Shield_07-Apr-2025_13-01-18';
fname_in = 'Shield_07-Apr-2025_13-01-18';

varnames = {'parsdata'}; % varaiables to be imported

% load aggregate data
for i = 1 : numel(varnames)
    load(fullfile(fdir_in, strcat(fname_in, '.mat')), varnames{i});
end

ii1 = [1, 2, 3]; % snapshots of post-flame agglomeration to be plotted

% baselines for locations of x-y-z axes in the plot of rendered aggregates
x_triad1 = [0.05, 0.46];
y0_triad1 = [0.13, 0.98];

% baselines for positions of tile titles (i.e. aggregate strcutural...
    % ...information)
x_ttl1 = [0.2, 0.6];
y0_ttl1 = [0.06, 0.96];

%% render selected aggregates and colorcode their primary particles...
    % ...based on shielding factor %%

% identify indices of aggregates to be plotted
jj1 = UTILS.MATCHAGG(vertcat(parsdata(ii1)), 'npp', 'sigmapp', 'n_hyb', ...
    0.8, 80, 50.0, 10000, 38, 300);  % x_max_raw = 300

n_agg_f1 = length(ii1); % number of selected snapshots

% display selected properties
for i = 1 : n_agg_f1
    fprintf('Group %d: npp = %d, sigmapp = %.4f, n_hyb = %d\n', i, ...
        parsdata(ii1(i)).npp(jj1(i)), ...
        parsdata(ii1(i)).sigmapp(jj1(i)), ...
        parsdata(ii1(i)).n_hyb(jj1(i)));
end

% initialize figure
f1 = figure(1);
f1.Position = [50, 50, 600, 900];
set(f1, 'color', 'white')
tl1_tot = tiledlayout(n_agg_f1,2);
tl1_tot.TileSpacing = 'compact';
tl1_tot.Padding = 'compact';

% generate colormap
opts1.cc = 'on';
opts1.cm = colormap("abyss");
opts1.cm = flip(opts1.cm,1);

% placeholders for tiles
tl1 = cell(n_agg_f1, 2);
row_ttl1 = cell(3,1);

% find triad position and title position
y_triad1 = linspace(y0_triad1(1), y0_triad1(2), n_agg_f1 + 1);
y_ttl1 = linspace(y0_ttl1(1), y0_ttl1(2), n_agg_f1 + 1);

% render aggregates
% for i = 1 : n_agg_f1
% 
%     % isometric view
%     tl1{i,1} = nexttile(tl1_tot, 2*i-1);
%     UTILS.PLOTPP_CONTINUOUS(parsdata(ii1(i)).pp{jj1(i)}(:,3),...
%         parsdata(ii1(i)).pp{jj1(i)}(:,4),...
%         parsdata(ii1(i)).pp{jj1(i)}(:,5),...
%         parsdata(ii1(i)).pp{jj1(i)}(:,2),...
%         parsdata(ii1(i)).spp{jj1(i)}, opts1);
%     clim([0 1])
%     hold on
%     UTILS.MAKEAX(f1, [x_triad1(1), y_triad1(i), 0.08, 0.08], '3d'); % triad
% 
%     % x-y plane
%     tl1{i,2} = nexttile(tl1_tot, 2*i);
%     UTILS.PLOTPP_CONTINUOUS(parsdata(ii1(i)).pp{jj1(i)}(:,3),...
%         parsdata(ii1(i)).pp{jj1(i)}(:,4),...
%         parsdata(ii1(i)).pp{jj1(i)}(:,5),...
%         parsdata(ii1(i)).pp{jj1(i)}(:,2),...
%         parsdata(ii1(i)).spp{jj1(i)}, opts1);
%     clim([0 1])
%     hold on
%     view(2)
%     UTILS.MAKEAX(f1, [x_triad1(2), y_triad1(i), 0.08, 0.08], 'xy'); % triad
% 
%     % generate title text 
%     row_ttl1{i} = strcat('$n_\mathrm{pp}$ =', {' '},...
%         num2str(parsdata(ii1(i)).npp(jj1(i))),...
%         ', $\sigma_\mathrm{pp}$ =', {' '},...
%         num2str(parsdata(ii1(i)).sigmapp(jj1(i)), '%.2f'),...
%         ', $n_\mathrm{hyb}$ =', {' '},...
%         num2str(parsdata(ii1(i)).n_hyb(jj1(i))));
% 
%     % print aggregate structural information
%     annotation('textbox', [x_ttl1(1), y_ttl1(i+1), x_ttl1(2), 0.03],...
%         'String', row_ttl1{i}, 'HorizontalAlignment', 'center',...
%         'VerticalAlignment', 'bottom', 'FontSize', 14, 'EdgeColor',...
%         'none', 'Interpreter', 'latex');
% 
% end
% 
% % generate colorbar showing shielding values
% cb1 = colorbar(tl1{1,1}, 'eastoutside');
% cb1.Layout.Tile = 'south';
% cb1.Label.String = 'Screening ratio [-]';
% cb1.FontSize = 12;
% cb1.Label.FontSize = 18;
% cb1.TickLabelInterpreter = 'latex';
% cb1.Label.Interpreter = 'latex';
% cb1.LineWidth = 1;

%% probability density functions of shielding factor over the course...
    % ...of post-flame agglomeration %%

% assign colors for pots-flame snapshots
clr2 = colormap(hot);
cind2 = round(1 + (length(clr2) - 1) .* (0.05 : 0.7 / (5 - 1) : 0.75)');
clr2 = clr2(cind2,:);
clr2(end,:) = [236,230,61] / 255;

% initialize figure
f2 = figure;
f2.Position = [100, 100, 500, 500];
set(f2, 'color', 'white');

n_shot = length(parsdata); % number of snapshots to be plotted

spp = cell(n_shot, 1); % allocate variable for ensemble shielding factors
mu_spp = zeros(n_shot, 1); % allocate means of shielding factor

% allocate probability density function (PDF) and evaluation points for...
    % ...ensemble shielding
f_spp = cell(n_shot, 1);
xi_spp = cell(n_shot, 1);

scale_spp = -0.25; % assign a scale for PFD to adjust curve extents

% allocate labels for x axis
xlbl2 = cell(n_shot, 1);

npp_tot = zeros(n_shot, 1); % allocate variable for the ensemble number...
    % ...of primary particles

for i = 1 : n_shot

    % concatinate shielding factors across aggregates
    spp{i} = cat(1, parsdata(i).spp{:});    
    
    % compute mean of shielding factor for each snapshot
    mu_spp(i) = mean(spp{i});

    % make a label for the distribution based on post-flame lifetime
    if i == 1
        xlbl2{i} = num2str(parsdata(i).r_n_agg(1), '%.0f');
    elseif ismember(i, [2,3])
        xlbl2{i} = num2str(parsdata(i).r_n_agg(1), '%.1f');
    else
        xlbl2{i} = num2str(parsdata(i).r_n_agg(1), '%.2f');
    end
    
    scatter(i-0.5, 1.2) % this is only to assign xticklabels
    hold on
    
    % print mean of each distribution
    text((i - 0.5), 1.03, ...
        sprintf('$\\overline{s}_\\mathrm{pp}$ = %.2f', mu_spp(i)), ...
        'Interpreter', 'latex', 'HorizontalAlignment', 'center',...       
        'FontSize', 10)

    % total number of primary particles across all aggregates
    npp_tot(i) = length(spp{i});
    
    % evaluate a probability density function od shielding factor for...
        % ...the entire population of primary particles
    [f_spp{i}, xi_spp{i}] = ksdensity(spp{i});
    
    % generate proper data format for shading beneath the distribution
    y_fill = [scale_spp * f_spp{i} + (i - 0.5),...
        (i - 0.5) * ones(size(f_spp{i}))];
    x_fill = [xi_spp{i}, fliplr(xi_spp{i})];
    
    plot(scale_spp * f_spp{i} + (i - 0.5), xi_spp{i}, 'Color', clr2(i,:),...
        'LineWidth', 2);
    fill(y_fill, x_fill, clr2(i,:), ...
        'FaceAlpha', 0.3, ...
        'EdgeColor', 'none')

end

xticks((1 : n_shot) - 0.5)  % specify tick positions for horizontal axis
xticklabels(xlbl2)  % assign labels to ticks

% set plot appearances
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02])
xlim([-0.3 4.8])
ylim([0 1])
xlabel('$n_\mathrm{agg}/(n_\mathrm{agg})_2$ [-]', 'interpreter', 'latex',...
    'FontSize', 18)
ylabel('$s_\mathrm{pp}$ [-]', 'interpreter', 'latex', 'FontSize', 18)

