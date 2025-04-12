clc
clear
clf('reset')
close all
warning('off')

%% initialize the script %%

% location of previously saved aggregate data (this should include...
    % ...spp field in pars strcuture for primary particle shielding)
fdir_in = 'F:\DLCA2\outputs\Shield_11-Apr-2025_00-35-50';
fname_in = 'Shield_11-Apr-2025_00-35-50';

varnames = {'parsdata'}; % varaiables to be imported

% load aggregate data
for i = 1 : numel(varnames)
    load(fullfile(fdir_in, strcat(fname_in, '.mat')), varnames{i});
end

ii1 = [1, 2, 3]; % snapshots of post-flame agglomeration to be plotted

% baselines for locations of x-y-z axes in the plot of rendered aggregates
x_triad1 = [0.07, 0.46];
y0_triad1 = [0.13, 0.98];

% baselines for positions of tile titles (i.e. aggregate strcutural...
    % ...information)
x_ttl1 = [0.2, 0.6];
y0_ttl1 = [0.06, 0.96];

% criteria on shielding ratio above which primary particles are...
    % ...not accounted when looking at projections of aggregates
spp_star = [1, 0.75, 0.5];

i_spp_star = 3; % selected criteria for plotting primary particle size...
    % ...vs aggregate size

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
opts.clim = [0 1];  % enforce full range from 0 to 1
opts.ft = 0.9;

% placeholders for tiles
tl1 = cell(n_agg_f1, 2);
row_ttl1 = cell(3,1);

% find triad position and title position
y_triad1 = linspace(y0_triad1(1), y0_triad1(2), n_agg_f1 + 1);
y_ttl1 = flip(linspace(y0_ttl1(1), y0_ttl1(2), n_agg_f1 + 1), 2);

% render aggregates
for i = 1 : n_agg_f1

    % isometric view
    tl1{i,1} = nexttile(tl1_tot, 2*i-1);
    UTILS.PLOTPP_CONTINUOUS(parsdata(ii1(i)).pp{jj1(i)}(:,3),...
        parsdata(ii1(i)).pp{jj1(i)}(:,4),...
        parsdata(ii1(i)).pp{jj1(i)}(:,5),...
        parsdata(ii1(i)).pp{jj1(i)}(:,2),...
        parsdata(ii1(i)).spp{jj1(i)}, opts1);
    % clim([0 1])
    hold on
    UTILS.MAKEAX(f1, [x_triad1(1), y_triad1(i), 0.08, 0.08], '3d'); % triad

    % x-y plane
    tl1{i,2} = nexttile(tl1_tot, 2*i);
    UTILS.PLOTPP_CONTINUOUS(parsdata(ii1(i)).pp{jj1(i)}(:,3),...
        parsdata(ii1(i)).pp{jj1(i)}(:,4),...
        parsdata(ii1(i)).pp{jj1(i)}(:,5),...
        parsdata(ii1(i)).pp{jj1(i)}(:,2),...
        parsdata(ii1(i)).spp{jj1(i)}, opts1);
    % clim([0 1])
    hold on
    view(2)
    UTILS.MAKEAX(f1, [x_triad1(2), y_triad1(i), 0.08, 0.08], 'xy'); % triad

    % generate title text 
    row_ttl1{i} = strcat('$n_\mathrm{pp}$ =', {' '},...
        num2str(parsdata(ii1(i)).npp(jj1(i))),...
        ', $\sigma_\mathrm{pp}$ =', {' '},...
        num2str(parsdata(ii1(i)).sigmapp(jj1(i)), '%.2f'),...
        ', $n_\mathrm{hyb}$ =', {' '},...
        num2str(parsdata(ii1(i)).n_hyb(jj1(i))));

    % print aggregate structural information
    annotation('textbox', [x_ttl1(1), y_ttl1(i), x_ttl1(2), 0.03],...
        'String', row_ttl1{i}, 'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'bottom', 'FontSize', 14, 'EdgeColor',...
        'none', 'Interpreter', 'latex');

end

% generate colorbar showing shielding values
cb1 = colorbar(tl1{1,1}, 'eastoutside');
cb1.Layout.Tile = 'south';
cb1.Label.String = '$S_\mathrm{pp}$ [-]';
cb1.FontSize = 12;
cb1.Label.FontSize = 18;
cb1.TickLabelInterpreter = 'latex';
cb1.Label.Interpreter = 'latex';
cb1.LineWidth = 1;

%% probability density functions of shielding factor over the course...
    % ...of post-flame agglomeration %%

% initialize figure
f2 = figure;
f2.Position = [100, 100, 500, 500];
set(f2, 'color', 'white');

% assign colors for pots-flame snapshots
clr2 = colormap(hot);
cind2 = round(1 + (length(clr2) - 1) .* (0.05 : 0.7 / (5 - 1) : 0.75)');
clr2 = clr2(cind2,:);
clr2(end,:) = [236,230,61] / 255;

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

hold on
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
    
    if i == 1
        % print mean of each distribution
        text((i - 0.85), 1.03, ...
            sprintf('$\\overline{S}_\\mathrm{pp}$ = %.2f', mu_spp(i)),...
            'Interpreter', 'latex', 'HorizontalAlignment', 'center',...       
            'FontSize', 12)
    else
        % print mean of each distribution
        text((i - 0.5), 1.03, num2str(mu_spp(i), '%.2f'),...
            'Interpreter', 'latex', 'HorizontalAlignment', 'center',...       
            'FontSize', 12)    
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
xlim([-0.5 5])
ylim([0 1])
xlabel('$n_\mathrm{agg}/(n_\mathrm{agg})_2$ [-]', 'interpreter', 'latex',...
    'FontSize', 18)
ylabel('$S_\mathrm{pp}$ [-]', 'interpreter', 'latex', 'FontSize', 18)

%% plot shielding factor vs. number of primary particles,...
    % ...colorcode based on polydispersity, and assign markers based...
    % ...on number of hybridity regions

% initialize figure
f3 = figure;
f3.Position = [150, 150, 500, 500];
set(f3, 'color', 'white');

% allocate number of aggregates in each snaphot
nagg = zeros(n_shot,1);

ms3 = [8, 16, 12, 16, 8]; % marker size
mt3 = {'^', 's', 'p', '*', 'o'}; % marker type

legtxt3 = cell(n_shot, 1); % allocate legends for post-flame snapshots

scat3 = cell(n_shot, 1); % allocate scatterplots

hold on
for i = 1 : n_shot
    
    % generate label for snapshots
    if i == 1
        legtxt3(i) = strcat('$n_\mathrm{agg}/(n_\mathrm{agg})_2$ =',...
            {' '}, num2str(parsdata(i).r_n_agg(1), '%.0f'));
    elseif ismember(i, [2,3])
        legtxt3(i) = strcat('$n_\mathrm{agg}/(n_\mathrm{agg})_2$ =',...
            {' '}, num2str(parsdata(i).r_n_agg(1), '%.1f'));
    elseif ismember(i, [4,5])
        legtxt3(i) = strcat('$n_\mathrm{agg}/(n_\mathrm{agg})_2$ =',...
            {' '}, num2str(parsdata(i).r_n_agg(1), '%.2f'));
    end

    nagg(i) = length(parsdata(i).pp); % find number of aggregates

    % allocate and calculate mean of shielding factor within...
        % ...individual aggregates
    parsdata(i).sagg = zeros(nagg(i),1);
    for k = 1 : nagg(i)
        parsdata(i).sagg(k) = mean(parsdata(i).spp{k});
    end
    
    scat3{i} = scatter(parsdata(i).npp, parsdata(i).sagg, ms3(i),...
        clr2(i,:), mt3{i}, 'LineWidth', 1);

end

% find bounds for x and y axes
bounds_npp_f3 = [0.7 * min(cat(1, parsdata.npp)),...
    1.3 * max(cat(1, parsdata.npp))];
bounds_sagg_f3 = [0.95 * min(cat(1, parsdata.sagg)),...
    1.05 * max(cat(1, parsdata.sagg))];

% set plot appearances
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim(bounds_npp_f3)
ylim(bounds_sagg_f3)
xlabel('$n_\mathrm{pp}$ [-]', 'interpreter', 'latex', 'FontSize', 18)
ylabel('$S_\mathrm{agg}$ [-]', 'interpreter', 'latex',...
    'FontSize', 18)

% print legends
legend(cat(1, scat3{:}), legtxt3, 'interpreter', 'latex',...
    'FontSize', 14, 'Location', 'southeast');

%% plot distributions of observed primary particle diameters %%

% number of shielding criteria
n_spp_star = length(spp_star);

% initialize figure
f4 = figure;
f4.Position = [200, 75, 400 * n_spp_star, 800];
set(f4, 'color', 'white');

% initialize layout
tl4 = tiledlayout(2, n_spp_star);
tl4.TileSpacing = 'compact';
tl4.Padding = 'compact';

% allocate observed primary particle diameter in 2d projections given a...
    % ...certain spp_star ((:,:,1) -> ensemble primary particle diameters,...
    % ...; (:,:,2) -> geometric mean of primary particle diameter within...
    % ...individual aggregates)
dpp_2d = cell(n_shot, n_spp_star, 2);

% allocate ensemble number of observed primary particle diameters
n_dpp_2d = zeros(n_shot, n_spp_star);

% ensemble geometric mean of primary particle diameter in 3d 
dpp_ens = cat(1, parsdata(1).pp{:});
dpp_ens = geomean(dpp_ens(:,2));

% allocate boxplots
bp4 = cell(n_shot, n_spp_star, 2);

% allocate subplot titles
titxt4 = cell(1, n_spp_star);

% allocate limits on y axes
bounds_dpp_f4 = zeros(2, n_spp_star, 2);

% extension factors on y axis limits
c_ylim = [0.1, 0.05];

% allocate decision variable for whether primary particles are (not)...
    % ...screened
kk = cell(n_shot, n_spp_star);

for i = 1 : n_shot
        
    for j = 1 : n_spp_star
        
        % allocate for individual aggregates
        dpp_2d{i,j,2} = zeros(nagg(i), 1);
        
        % allocate shielding factor decision variable for each snapshot
        kk{i,j} = cell(nagg(i),1);

        % filter primary particles based on spp_star and calculate...
            % ...geometric means of their diameter
        for k = 1 : nagg(i)
            kk{i,j}{k} = parsdata(i).spp{k} <= spp_star(j);
            dpp_2d{i,j,2}(k) = geomean(parsdata(i).pp{k}(kk{i,j}{k},2));
        end
    
        % compile primary particles acrosss all aggregates
        dpp_2d{i,j,1} = cat(1, parsdata(i).pp{:});

        % filter ensemble of primary particles based on shielding...
            % ...factor and calculate geometric mean
        dpp_2d{i,j,1} = dpp_2d{i,j,1}(cat(1, kk{i,j}{:}),2);
        
        % total number of observable primary particles
        n_dpp_2d(i,j) = length(dpp_2d{i,j,1});
    
        for l = 1 : 2

            nexttile(j + (l-1) * n_spp_star)

            % plot distribution of primary particle diameter
            bp4{i,j,l} = boxplot(1e9 * dpp_2d{i,j,l}, 'Positions', i,...
                'Notch', 'on', 'Symbol', 'o', 'Widths', 0.5);
            hold on
        
            % Find the box object
            boxObj = findobj(bp4{i,j,l}, 'Tag', 'Box');    
            
            % fill inside the box
            patch(get(boxObj, 'XData'), get(boxObj, 'YData'), clr2(i, :),...
                'FaceAlpha', 0.3, 'EdgeColor', clr2(i, :), 'LineWidth', 1);
            
            % adjust the median line
            boxMed = findobj(bp4{i,j,l}, 'Tag', 'Median');
            set(boxMed, 'Color', clr2(i, :), 'LineWidth', 2);
            
            % adjust outlier markers
            outliers = findobj(bp4{i,j,l}, 'Tag', 'Outliers');
            outliers.MarkerEdgeColor = clr2(i, :);
            outliers.MarkerSize = 3;
            
            % adjust whiskers
            upwhisker = findobj(bp4{i,j,l},'type', 'line', 'tag', 'Upper Whisker');
            set(upwhisker, 'linestyle', '-');
            lowwhisker= findobj(bp4{i,j,l}, 'type', 'line','tag', 'Lower Whisker');
            set(lowwhisker, 'linestyle', '-');
    
            xticks(1 : n_shot)  % specify tick positions for horizontal axis
            xticklabels(xlbl2)  % assign labels to ticks
            
            if i == n_shot
    
                % % calculate and plot ensemble mean
                % kkk = cat(1, parsdata(i).spp{:}) <= spp_star(j);
                % dpp_ens(j) = geomean(pp_ens(kkk,2));
                plt4_ens = plot(linspace(0, n_shot+1, 100),...
                    1e9 * repmat(dpp_ens,1,100), 'Color', [0, 0, 0],...
                    'LineWidth', 1.5, 'LineStyle', ':');
            
                % set plot appearances
                box on
                set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
                    'TickLength', [0.02 0.02], 'YScale', 'log')
                xlim([0.25 n_shot+0.75])
                
                % get ymin and ymax
                bounds_dpp_f4(1,j,l) = min(1e9 * cat(1,dpp_2d{:,j,l}));
                bounds_dpp_f4(2,j,l) = max(1e9 * cat(1,dpp_2d{:,j,l}));
                
                % assign label for vertical axes
                if j == 1
                    if l ==1
                        ylabel('$d_\mathrm{pp}^{(i)}$ [nm]',...
                            'interpreter', 'latex', 'FontSize', 18)
                    else
                        ylabel('$d_\mathrm{pp}$ [nm]',...
                            'interpreter', 'latex', 'FontSize', 18)
                    end
                end
            end

            nexttile(j)

            % make a title for vertical tiles showing enforced spp_star
            titxt4{j} = sprintf('$S_\\mathrm{pp}^\\mathrm{*}$ = %.2f',...
                spp_star(j));
            title(titxt4{j}, 'interpreter', 'latex', 'FontSize', 14)

            % add some space behind title
            subtitle(' ', 'interpreter', 'latex', 'FontSize', 10)

            yticks(cat(2, linspace(1,10,10), linspace(20,100,9)))

        end
        
    end

end

xlabel(tl4, '$n_\mathrm{agg}/(n_\mathrm{agg})_2$ [-]', 'interpreter',...
    'latex', 'FontSize', 18) % label for horizontal axis

% legend showing ensemble diameter for the lower tiles
legend(plt4_ens, '$\langle{d_\mathrm{pp}}\rangle$',...
    'interpreter', 'latex', 'FontSize', 14, 'location', 'southeast')

% set identical limits on y axis
for j = 1 : n_spp_star
    for l = 1 : 2
        nexttile((l-1) * n_spp_star + j)      
        ylim([(1 - c_ylim(l)) * min(bounds_dpp_f4(1,:,l)),...
            (1 + c_ylim(l)) * max(bounds_dpp_f4(2,:,l));])
    end
end

%% plot mean primary particle diameter vs. aggregate area diameter...
    % ..."considering shielding" %%

% initialize figure
f5 = figure(5);
f5.Position = [250, 125, 800, 500];
set(f5, 'color', 'white')

% initialize layout
tl5 = tiledlayout(1, 2);
tl5.TileSpacing = 'compact';
tl5.Padding = 'compact';

plt5 = cell(n_shot + 2, 1); % initialize placholders for plots

% assign variables for universal correlation
D_TEM = 0.35; % exponent
dpp_100 = 17.8; % pefactor
da_lim_uc = [1e0 2e4];  % limits on the projected area diameter
n_da_uc = 1e4; % number of data
uc1 = @(y) dpp_100 * (y / 100) .^ D_TEM; % on-demand function for the...
    % ...forward correlation in the geometrical domain (dpp as a...
    % ...function of da in [nm])
r_uc1 = (da_lim_uc(2) / da_lim_uc(1)) ^ (1 / (n_da_uc - 1));
da_uc = da_lim_uc(1) * ones(n_da_uc,1) .* r_uc1 .^ (((1 : n_da_uc) - 1)');
dpp_uc = uc1(da_uc);

nexttile(1)

% plot universal correlation of of Olfert & Rogak (2019)
plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

% plot ensemble primary particle diameter
plot(da_uc, 1e9 * repmat(dpp_ens, size(da_uc)),...
    'Color', [0, 0, 0], 'LineStyle', ':', 'LineWidth', 2);

% plot dpp (2D projections) vs da over post-flame snapshots
for i = 1 : n_shot
    scatter(1e9 * parsdata(i).da, 1e9 * dpp_2d{i,1,2},...
        ms3(i), clr2(i,:), mt3{i}, 'LineWidth', 1);
end

% set plot appearances
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 18)
ylabel('$d_\mathrm{pp}^\mathrm{(2D)}$ [nm]', 'interpreter',...
    'latex', 'FontSize', 18)

% apply proper bounds to x and y axes
bounds_da_f5 = [1e9 * 0.8 * min(cat(1, parsdata.da)),...
    1e9 * 1.2 * max(cat(1, parsdata.da))];
bounds_dpp_f5 = [1e9 * 0.95 * min(cat(1, dpp_2d{:,3,2})),...
    1e9 * 1.05 * max(cat(1, dpp_2d{:,3,2}))];xlim(bounds_da_f5)
ylim(bounds_dpp_f5)

% title stating criteria for shielding factor
title(sprintf('$S_\\mathrm{pp}^\\mathrm{*}$ = %.2f', spp_star(1)),...
    'interpreter', 'latex', 'FontSize', 14)
% add some space behind title
subtitle(' ', 'interpreter', 'latex', 'FontSize', 8)

nexttile(2)

% plot universal correlation of of Olfert & Rogak (2019)
plt5{end} = plot(da_uc, dpp_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 3);
hold on

% plot ensemble primary particle diameter
plt5{end-1} = plot(da_uc, 1e9 * repmat(dpp_ens, size(da_uc)),...
    'Color', [0, 0, 0], 'LineStyle', ':', 'LineWidth', 2);

% plot dpp (2D projections) vs da over post-flame snapshots
for i = 1 : n_shot
    plt5{i} = scatter(1e9 * parsdata(i).da, 1e9 * dpp_2d{i,i_spp_star,2},...
        ms3(i), clr2(i,:), mt3{i}, 'LineWidth', 1);
end

% set plot appearances
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 18)

% apply proper bounds to x and y axes
xlim(bounds_da_f5)
ylim(bounds_dpp_f5)

% title stating criteria for shielding factor
title(sprintf('$S_\\mathrm{pp}^\\mathrm{*}$ = %.2f', spp_star(i_spp_star)),...
    'interpreter', 'latex', 'FontSize', 14)
% add some space behind title
subtitle(' ', 'interpreter', 'latex', 'FontSize', 8)

lgd5 = legend(cat(1, plt5{:}), cat(1,legtxt3,...
    {'$\langle{d}_\mathrm{pp}\rangle$'}, {'Olfert $\&$ Rogak (2019)'}),...
    'interpreter', 'latex', 'FontSize', 14);
lgd5.Layout.Tile = 'south';

%% save plots and workspace %%

% make a directory to save outputs
dir0_out = datestr(datetime('now'));
dir0_out = regexprep(dir0_out, ':', '-');
dir0_out = regexprep(dir0_out, ' ', '_');
dir_out = strcat('outputs\', 'PostShield_', dir0_out, '\');
if ~isfolder(dir_out)
    mkdir(dir_out); % if it doesn't exist, create the directory
end

% save worksapce
save(strcat(dir_out, 'PostShield_', dir0_out, '.mat'))

% print figures
exportgraphics(f1, strcat(dir_out, 'render-shield.jpg'),...
    'BackgroundColor','none', 'Resolution', 300)
exportgraphics(f2, strcat(dir_out, 'dist-shield-pp.jpg'),...
    'BackgroundColor','none', 'Resolution', 300)
exportgraphics(f3, strcat(dir_out, 'scat-shield-agg.jpg'),...
    'BackgroundColor','none', 'Resolution', 300)
exportgraphics(f4, strcat(dir_out, 'dpp-bias.jpg'),...
    'BackgroundColor','none', 'Resolution', 300)
exportgraphics(f5, strcat(dir_out, 'dpp-vs-da-bias.jpg'),...
    'BackgroundColor','none', 'Resolution', 300)
