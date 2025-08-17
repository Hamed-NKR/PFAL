function h = PCF_MULTI_v2(g_multi, r_multi, pp_multi, opts)
% PCF_MULTI gets multiple pair-correlation function (PCF) distributions...
%   ...from different scenarios and organizes them for comparison...
%   ...in a table of PCF vs. radial location plots. The plots...
%   ...are separated based on the number of primaries, polydispersity...
%   ...and hybridity of aggregates/agglomerates.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   g_multi: An m*n*l cell array containing PCF data. Each cell contains...
%       ...data from a certain trial of a certain scenario of hybridity;...
%       ...i.e, (m) denotes hybridity level, (n) denotes the controlled...
%       ...polydispersity and size, and (l) denotes whether we look at...
%       ...data from a real hybrid agglomerate or one of the repeats of...
%       ...a uniform aggregate with similarly initialized parameters...
%       ...(number of primaries, geometric standard deviation of primary
%       ...particle size).
%   r_multi: An m*n*l cell array including the raidal logarithmic points...
%       ...of calculation corresponding to the data in g_multi.
%   pp_multi: An m*n*l cell array of primary particle info matrices...
%       ...corresponding to the g_multi+r_multi data.
%   opts: visibility options
% ----------------------------------------------------------------------- %
% 
% Outputs:
%   h: PCF plot figure handle
% ----------------------------------------------------------------------- %

% initialize the structure for visibility feature if not existing 
if ~exist('opts', 'var') 
    opts = struct();
end

if ~isfield(opts, 'sz')
    opts.sz = []; % initialize field for the size of subplots
end
if isempty(opts.sz)
    opts.sz = 500; % set default size for subplots
end

% number of rows and columns in subplots and number of lines within them
[m,n,l] = size(g_multi);

if ~isfield(opts, 'dimx')
    opts.dimx = []; % initialize field for horizontal location of textboxes
end
if isempty(opts.dimx)
    % set default horizontal box locations
    if n == 3
        opts.dimx = [0.072 0.380 0.687];
    end
end

if ~isfield(opts, 'dimy')
    opts.dimy = []; % initialize field for vertical location of textboxes
end
if isempty(opts.dimy)
    % set default vertical box locations
    if m == 2
        opts.dimy = [0.413 0];
    end
end

if ~isfield(opts, 'dimyy')
    opts.dimyy = []; % initialize field for the secondary vertical adjustment
end
if isempty(opts.dimyy)
    % set default for vertical adjustment value
    if m == 2
        opts.dimyy = 0.061;
    end
end
% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, n * opts.sz, m * opts.sz];
set(h, 'color', 'white');

% initialize layout
tt = tiledlayout(m,n);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

lin_c1 = [203,24,29;
    33,113,181;
    35,139,69;
    106,81,163;
    217,71,1;
    82,82,82] / 255; % line color repository - stream 1

lin_c2 = [252,174,145;
    189,215,231;
    186,228,179;
    203,201,226;
    253,190,133;
    204,204,204] / 255; % line color repository - stream 2

lin_w = [3, 1.5]; % line width repo.
lin_t = {'-', ':'}; % line style repo.

% Placeholder for aggregate mass and gyration diameter
dg = zeros(10,1);
ma = zeros(10,1);

for i = 1 : m
    for j = 1 : n
        nexttile
        
        % plot the hybrid agglomerates
        plt1 = plot(r_multi{i,j,1}, g_multi{i,j,1}, 'Color', lin_c1(i,:),...
            'LineStyle', lin_t{1}, 'LineWidth', lin_w(1));
        
        hold on
        
        % make the subtitles with properties of agglomerates
        nhyb = length(unique(pp_multi{i,j,1}(:,6)));
        sigma = UTILS.GEOSTD(pp_multi{i,j,1}(:,2));
        npp = size(pp_multi{i,j,1},1);
        ttl = strcat({'$n_\mathrm{hyb}$ ='}, {' '}, num2str(nhyb, '%d'),...
            {','}, ' $\sigma_\mathrm{pp}$ =', {' '}, num2str(sigma, '%.3f'),...
            {','}, {' $n_\mathrm{pp}$ ='}, {' '}, num2str(npp, '%d'));
        
        % calculate mass, gyration diameter and mean pp size to ensure case control
        dg1 = PAR.GYRATION(pp_multi(i,j,1),npp);
        ma1 = (1860 * pi / 6) * sum(pp_multi{i,j,1}(:,2).^3);
        dpp1 = geomean(pp_multi{i,j,1}(:,2));
        
        % plot and analyze the equivalent uniform aggregates
        for k = 2 : l
            plt2 = plot(r_multi{i,j,k}, g_multi{i,j,k}, 'Color', lin_c2(i,:),...
                'LineStyle', lin_t{2}, 'LineWidth', lin_w(2));
            
            % normalize the size of primaries in uniform aggeragtes
            dpp2 = geomean(pp_multi{i,j,k}(:,2));
            rpp = dpp1 / dpp2;
            pp_multi{i,j,k}(:,2) = rpp * pp_multi{i,j,k}(:,2);
            
            dg(k-1) = PAR.GYRATION(pp_multi(i,j,k),npp);
            ma(k-1) = (1860 * pi / 6) * sum(pp_multi{i,j,k}(:,2).^3);
        end
                
        % label hybrids and non-hybrids
        if (j == n)
            legend([plt1, plt2], {'Hybrid', 'Uniform'},...
                'interpreter', 'latex', 'FontSize', 18,...
                'Location', 'northeast');
        end
        
        % average aggregate mass and gyration diameter for the uniform set
        dg2 = mean(dg);
        ci_dg2 = 1.96 * std(dg) / sqrt(l-1);
        ma2 = mean(ma);
        ci_ma2 = 1.96 * std(ma) / sqrt(l-1);
        
        % make the annotation box
        dim = [opts.dimx(j), opts.dimy(i), 0.3 , 0.3];
        str = strcat('for', {' '}, '$d_\mathrm{pp}$ =', {' '},...
                num2str(1e9 * dpp1, '%.1f'), ' nm:',...
            string(newline), '$\left(d_\mathrm{g}\right)_\mathrm{hyb}$ =', {' '},...
                num2str(1e9 * dg1, '%.1f'), ' (nm)',...
            string(newline), '$\left(\overline{d}_\mathrm{g}\right)_\mathrm{uni}$ =',...
                {' '}, num2str(1e9 * dg2, '%.1f'), {'   '}, {'$\pm$'},...
                {' '}, num2str(1e9 * ci_dg2, '%.1f'), ' (nm)',...
            string(newline), '$\left(m_\mathrm{agg}\right)_\mathrm{hyb}$ =', {' '},...
                num2str(ma1, '%.2e'), ' (kg)',...
            string(newline), '$\left(\overline{m}_\mathrm{agg}\right)_\mathrm{uni}$ =',...
                {' '}, num2str(ma2, '%.2e'), {' '}, {'$\pm$'},...
                string(newline), {'   '}, num2str(ci_ma2, '%.2e'), ' (kg)');
        ant = annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on',...
            'interpreter', 'latex', 'FontSize', 14);
        if i > 1; p = ant.Parent; p.Children(1).Position(2) = ...
                p.Children(1).Position(2) - opts.dimyy(i-1); end
        
        box on
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18,...
            'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
        xlim([0.1 130])
        ylim([3e-5 1])
        ax = gca;
        set(ax, 'xtick',[0.1, 1, 10, 100])
        set(ax, 'ytick',[1e-4, 1e-3, 1e-2, 1e-1, 1])
        if ismember(j, 2:n)
            set(ax, 'yticklabel',[])
        end
        if ismember(i, 1:(m-1))
            set(ax, 'xticklabel',[])
        end
        
        title(ttl, 'FontSize', 20, 'interpreter','latex')
    end
end

xlabel(tt, '$\hat{r}$ (-)', 'interpreter', 'latex', 'FontSize', 20)
ylabel(tt, '$g\left(\hat{r}\right)$ (-)', 'interpreter', 'latex', 'FontSize', 20)

end

