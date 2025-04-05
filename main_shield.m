% ---------------- MAIN SCRIPT ---------------- %
clc
clear
clf('reset')
close all
warning('off')

%% initialize the script %%
fdir_in = 'F:\DLCA2\outputs\postLD2-01-Apr-2025_11-54-39_LD2-25NOV24';
fname_in = 'Post_LD2-25NOV24';

varnames = {'parsdata'};
for i = 1 : numel(varnames)
    load(fullfile(fdir_in, strcat(fname_in, '.mat')), varnames{i});
end

ii0 = [1, 2, 3];
n_agg_plt = length(ii0);

triad_pos = {[0.07, 0.71, 0.08, 0.08], [0.49, 0.71, 0.08, 0.08];
             [0.07, 0.42, 0.08, 0.08], [0.49, 0.42, 0.08, 0.08];
             [0.07, 0.11, 0.08, 0.08], [0.49, 0.11, 0.08, 0.08]};
y_row_ttl = [0.96, 0.66, 0.37];
x_row_ttl = [0.2, 0.6];

% call MATCHAGG with fixed seed to ensure reproducibility and cap on npp
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

opts2.cc = 'on';
opts2.cm = colormap("abyss");
opts2.cm = flip(opts2.cm,1);

tl2 = cell(n_agg_plt, 2);
row_ttl = cell(3,1);

for i = 1 : n_agg_plt
    tl2{i,1} = nexttile(tl2_tot, 2*i-1);
    UTILS.PLOTPP(parsdata(ii0(i)).pp{jj0(i)}(:,3), parsdata(ii0(i)).pp{jj0(i)}(:,4), ...
                 parsdata(ii0(i)).pp{jj0(i)}(:,5), parsdata(ii0(i)).pp{jj0(i)}(:,2), ...
                 parsdata(ii0(i)).pp{jj0(i)}(:,2), opts2);
    hold on
    UTILS.MAKEAX(f2, triad_pos{i,1}, '3d');

    tl2{i,2} = nexttile(tl2_tot, 2*i);
    UTILS.PLOTPP(parsdata(ii0(i)).pp{jj0(i)}(:,3), parsdata(ii0(i)).pp{jj0(i)}(:,4), ...
                 parsdata(ii0(i)).pp{jj0(i)}(:,5), parsdata(ii0(i)).pp{jj0(i)}(:,2), ...
                 parsdata(ii0(i)).pp{jj0(i)}(:,2), opts2);
    hold on
    view(2)
    UTILS.MAKEAX(f2, triad_pos{i,2}, 'xy');

    row_ttl{i} = strcat('$n_\mathrm{pp}$ =', {' '}, num2str(parsdata(ii0(i)).npp(jj0(i))), ...
        ', $\sigma_\mathrm{pp}$ =', {' '}, num2str(parsdata(ii0(i)).sigmapp(jj0(i)), '%.2f'), ...
        ', $n_\mathrm{hyb}$ =', {' '}, num2str(parsdata(ii0(i)).n_hyb(jj0(i))));

    annotation('textbox', [x_row_ttl(1), y_row_ttl(i), x_row_ttl(2), 0.03], ...
        'String', row_ttl{i}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 14, 'EdgeColor', 'none', 'Interpreter', 'latex');
end

cb2 = colorbar(tl2{1,1}, 'eastoutside');
cb2.Layout.Tile = 'south';
cb2.Label.String = 'Shielding ratio [-]';
cb2.FontSize = 12;
cb2.Label.FontSize = 18;
cb2.TickLabelInterpreter = 'latex';
cb2.Label.Interpreter = 'latex';
cb2.LineWidth = 1;
