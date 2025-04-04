clc
clear
clf('reset')
close all
warning('off')

%% initialize the script %%

% address of second-stage langevin dynamics data to be imported
fdir_in = 'F:\DLCA2\outputs\postLD2-01-Apr-2025_11-54-39_LD2-25NOV24';
fname_in = 'Post_LD2-25NOV24';

% variables of interest in the data
varnames = {'parsdata'};

% load second-stage LD aggregates
for i = 1 : numel(varnames)
    load(strcat(fdir_in, '\', fname_in, '.mat'), varnames{i})
end

% assign a subset of aggregates for plotting

ii0 = [1, 3, 5]; % determine post-flame snapshots from which the...
    % ...aggregates to be plotted are chosen

% set target values on number of primary particles and polydispersity...
    % ...to choose aggregates being plotted
nPP_target = 110;
sigmaPP_target = 1.3;

% fixed triad locations (one per tile)
triad_pos = {[0.08, 0.9, 0.08, 0.08], [0.85, 0.9, 0.08, 0.08];
    [0.08, 0.43, 0.08, 0.08], [0.55, 0.43, 0.08, 0.08];
    [0.08, 0.3, 0.08, 0.08], [0.85, 0.3, 0.08, 0.08]};
    
%% find aggregates to be rendered %%

n_agg_plt = length(ii0); % numebr of hybridity snapshots being rendered

% normalized difference in number of primary particles between...
    % ...repository and target
dnPP_hat = cell(n_agg_plt,1);

% normalized difference in polydispersity
dsigmaPP_hat = cell(n_agg_plt,1);

% indices of aggregates being rendered
jj0 = zeros(n_agg_plt,1);

% calculate normalized differences
for i = 1 : n_agg_plt

    dnPP_hat{i} = abs(parsdata(ii0(i)).npp - nPP_target); 
    dnPP_hat{i} = dnPP_hat{i} / max(dnPP_hat{i});

    dsigmaPP_hat{i} = abs(parsdata(ii0(i)).sigmapp - sigmaPP_target); 
    dsigmaPP_hat{i} = dsigmaPP_hat{i} / max(dsigmaPP_hat{i});
    
    % find aggregate that is closet to target (via minimizing differences)
    jj0(i) = find((dnPP_hat{i} + dsigmaPP_hat{i}) ==...
        min(dnPP_hat{i} + dsigmaPP_hat{i}), 1);

end

%% plot subselected aggregates %%

% initialize temporal dpp vs da figure
f2 = figure(2);
f2.Position = [100, 100, 600, 900];
set(f2, 'color', 'white')

% initialize layout
tl2_tot = tiledlayout(n_agg_plt,2);
tl2_tot.TileSpacing = 'compact';
tl2_tot.Padding = 'compact';

% define rendering color
opts2.cc = 'on';
opts2.cm = colormap("abyss");
opts2.cm = flip(opts2.cm,1);

tl2 = cell(n_agg_plt, 2); % placeholder for individual tiles

for i = 1 : n_agg_plt

    % first row: a nonhybrid aggregate
    % other rows: aggregates with increasing hybridity (as moving down...
        % ...the rows) but with the same number of primary particles and...
        % ...polydispersity

    tl2{i,1} = nexttile(tl2_tot, 2*i-1); % isometric view
    
    UTILS.PLOTPP(parsdata(ii0(1)).pp{jj0(i)}(:,3), parsdata(ii0(1)).pp{jj0(i)}(:,4),...
        parsdata(ii0(1)).pp{jj0(i)}(:,5), parsdata(ii0(1)).pp{jj0(i)}(:,2),...
        parsdata(ii0(1)).pp{jj0(i)}(:,2), opts2); % plot aggregate
    hold on
    
    UTILS.MAKEAX(f2, triad_pos{i,1}, '3d'); % render x-y-z axes
    
    tl2{i,2} = nexttile(tl2_tot, 2*i); % x-y view
    
    UTILS.PLOTPP(parsdata(ii0(1)).pp{jj0(i)}(:,3), parsdata(ii0(1)).pp{jj0(i)}(:,4),...
        parsdata(ii0(1)).pp{jj0(i)}(:,5), parsdata(ii0(1)).pp{jj0(i)}(:,2),...
        parsdata(ii0(1)).pp{jj0(i)}(:,2), opts2); % plot aggregate
    hold on
    view(2)

    UTILS.MAKEAX(f2, triad_pos{i,2}, 'xy'); % render x-y-z axes

end

% add a shared colorbar under the whole layout
cb2 = colorbar(tl2{1,1}, 'eastoutside');  % generate colorbar in first tile
cb2.Layout.Tile = 'south'; % move it to the bottom
cb2.Label.String = 'Shielding ratio [-]'; % add label
cb2.FontSize = 12; % tick label font size
cb2.Label.FontSize = 18; % label font size
cb2.TickLabelInterpreter = 'latex'; % tick label format
cb2.Label.Interpreter = 'latex'; % label format
cb2.LineWidth = 1; % width of edges and ticks



    