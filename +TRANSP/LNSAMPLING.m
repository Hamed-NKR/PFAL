function [d, ind, fn, fn0, h] = LNSAMPLING(d0, mu_d, sigma_d, n_agg,...
    n_bin, del_d, opts)
    
% "LNSAMPLING" samples a lognormal size distribution from a population...
    % ...of aggregates of different sizes.
% ----------------------------------------------------------------------- %
%
% Inputs:
%   d0: Pre-sampling diameter set of a given aggregate population
%   mu_d: The desired mean size after sampling
%   sigma_d: The desired sampling standard deviation
%   n_agg: Number of target particles to be selected
%   n_bin: Number of bins for discretization
%   del_d: Range of sampling
%   opts: Function options
% ----------------------------------------------------------------------- %
% 
% Outputs:
%   d: Post-sampling lognormally distributed diameter set
%   inds: Indices of particles selected
%   fn0: Number frequency of pre-sampling size dataset with respect to...
%       ...the bins defined
%   fn: Number frequency of post-sampling
%   h: Figure handle for size distribution plots
% ----------------------------------------------------------------------- %

n_agg0 = length(d0); % initial number of aggregates

% get the size range
if ~exist('del_d', 'var') || isempty(del_d)
    [del_d(1), del_d(2)] = bounds(d0);
else
    del_d = sort(del_d(1:2));
end

% set size bins
r_d = (del_d(2) / del_d(1))^(1 / n_bin);
d_bin = del_d(1) * ones(n_bin + 1, 1);
for i = 1 : n_bin
    d_bin(i+1) = d_bin(i+1) * r_d^(i);
end
d_c = zeros(n_bin,1); % bin centers

pl = zeros(n_agg0,1); % initialize particle bin label set (range: 1-n_bin)

% initialize number of output particles, size mean, std, number of bins,...
    % ...and visualization options if not given
if (~exist('n_agg', 'var'))
    n_agg = [];
end

if (~exist('mu_d', 'var')) || isempty(mu_d) || (~isnumeric(mu_d))
%     mu_d = exp(median(log(sort(d0)))); % default to be...
%         % ...logarithmic median
    mu_d = sqrt(del_d(1) * del_d(2)); % default to be...
        % ...logarithmic mean of max and mean
end

if (~exist('sigma_d', 'var')) || isempty(sigma_d) || (~isnumeric(sigma_d))
    sigma_d = 1.4;
end

if (~exist('n_bin', 'var')) || isempty(n_bin) || (~isnumeric(n_bin))
    n_bin = 10;
end

if ~exist('opts', 'var') 
    opts = struct();
end

% set visualization option
if ~isfield(opts, 'visual')
    opts.visual = [];
end
opts_visual = opts.visual;
if isempty(opts_visual)
    opts_visual = 'off'; % default not to plot the outputs
end

% set aggregate size type
if ~isfield(opts, 'randvar')
    opts.randvar = [];
end
opts_randvar = opts.randvar;

% set the data fitting type
if ~isfield(opts, 'nfit')
    opts.nfit = [];
end
opts_nfit = opts.nfit;
if isempty(opts_nfit)
    opts_nfit = 'off'; % default to fit based on the total number of aggs (not the lognormal peak)
end

% assign the xlabel for the output graphs
if ismember(opts_visual, {'ON', 'On', 'on'})
    x_lb = "d";
    switch opts_randvar
        case {'MOBILITY', 'Mobility', 'mobility'}
            x_lb = x_lb + "_m";
        case {'AREA', 'Area', 'area'}
            x_lb = x_lb + "_a";
    end
end

% initialize fequency arrays
fn0 = zeros(n_bin, 1);
fn = zeros(n_bin, 1);
fn_fit = zeros(n_bin, 1);

d = cell(n_bin,1); % initialize post-sampling size array
ind = cell(n_bin,1); % placeholder for the aggs selected by LN sampling

% put aggregates in their bins
for i = 1 : n_bin
    if i == 1 % The first bin (only upper limit)
        ii = d0 < d_bin(i+1);
    elseif i < n_bin
        ii = (d0 >= d_bin(i)) & (d0 < d_bin(i+1));
    else % the last bin (all the remaining)
        ii = pl == 0;
    end
    
    d_c(i) = sqrt(d_bin(i) * d_bin(i+1));
    pl(ii) = i; % assign bin labels to the aggregates
    fn0(i) = nnz(ii) / n_agg0; % size frequency of pre-sampling population   
    
    fn_fit(i) = (normcdf(log(d_bin(i+1)), log(mu_d), log(sigma_d)) -...
        normcdf(log(d_bin(i)), log(mu_d), log(sigma_d))); % the lognormal fit to the...
            % ...frequency data
    
    ind{i} = find(pl == i); % initialize i'th bin post-sampling indices
    d{i} = d0(pl == i); % ~ diameters
end

% find the number of output aggregates
if ismember(opts_nfit, {'ON', 'On', 'on'})
    jmax = find(d_bin > mu_d, 1) - 1;
    n_agg = ceil(fn0(jmax) * n_agg0 / fn_fit(jmax));
    if ~isempty(n_agg)
        warning('Number of aggregates adjusted for lognormal sampling...')
    end
elseif isempty(n_agg)
    n_agg = n_agg0;
end

% randomly select from the bins
for i = 1 : n_bin
    n_fit = ceil(fn_fit(i) * n_agg); % number of particles to be...
        % ...selected in each bin
    n0 = round(fn0(i) * n_agg0); % number of pre-sampling particles in each bin
    
    if n_fit <= n0 % if enough particles exist in bin for random sampling
        iii= randperm(n0,n_fit); % index of particles in bin to be sampled
        
    elseif n0 > 0
        c_rep = ceil(n_fit / n0); % data replication factor
        
        % replicate data to compensate for the missing frequency
        ind{i} = repmat(ind{i}, c_rep, 1);
        d{i} = repmat(d{i}, c_rep, 1);
        
        iii = randperm(c_rep * n0, n_fit); % adjust the bin indices for...
            % ...replication        
    end
    
    % assign the data of selected particles 
    if n0 ~= 0
        ind{i} = ind{i}(iii);
        d{i} = d{i}(iii);
        
        fn(i) = length(d{i}) / n_agg; % get post-sampling frequency
    end 
end

% plot the pre- and post-sampling distributions if requested
if ismember(opts_visual, {'ON', 'On', 'on'})
    % initialize figure properties
    figure
    h = gcf;
    if ~all(h.Position == [0, 0, 1200, 600])
        h.Position = [0, 0, 1200, 600];
    end
    set(h, 'color', 'white')
    
    % set the layout
    tt = tiledlayout(1, 2);
    tt.TileSpacing = 'compact';
    tt.Padding = 'compact';
    
    tt1 = nexttile;
    
    % initial size distribution
    x0 = 1e9 * [d_bin(1); repelem(d_bin(2 : end - 1), 2); d_bin(end)];
    dfn0_dlogda = fn0 ./ log(d_bin(2 : end) ./ d_bin(1 : end - 1));
    y11 = repelem(dfn0_dlogda,2);
    p11 = plot(x0, y11, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2.5);
    hold on
    
    % desired size distribution
    x12 = (1e9 * exp(linspace(log(del_d(1)), log(del_d(2)), 1000)))';
%     y12 = 100 * exp(-(log(x12) - log(1e9 * mu_d)).^2 / (2 * (log(sigma_d))^2)) ./...
%         (sqrt(2 * pi) * log(sigma_d) * log(x12));
    y12 = pdf('Normal', log(x12), log(1e9 * mu_d), log(sigma_d));
    p12 = plot(x12, y12, 'Color', [0.4660 0.6740 0.1880],...
        'LineWidth', 2.5, 'LineStyle', '-.');
    
    box on
    set(gca, 'FontName', 'SansSerif', 'FontSize', 12, 'TickLength', [0.02 0.02])
    set(gca, 'XScale', 'log')
    title(tt1, 'Before sampling', 'FontName', 'SansSerif',...
        'FontWeight', 'bold', 'FontSize', 16)
    
    tt2 = nexttile;
    
    % final size distribution
    dfn_dlogda = fn ./ log(d_bin(2 : end) ./ d_bin(1 : end - 1));
    y2 = repelem(dfn_dlogda,2);
    p21 = plot(x0, y2, 'Color', [0 0.4470 0.7410], 'LineWidth', 2.5);
    hold on
    
    plot(x12, y12, 'Color', [0.4660 0.6740 0.1880],...
        'LineWidth', 2.5, 'LineStyle', '-.')
    
    box on
    set(gca, 'FontName', 'SansSerif', 'FontSize', 12, 'TickLength', [0.02 0.02])
    set(gca, 'XScale', 'log')
    title(tt2, 'After sampling', 'FontName', 'SansSerif',...
        'FontWeight', 'bold', 'FontSize', 16)
    
    xlabel(tt, x_lb + ' (nm)', 'FontName', 'SansSerif', 'FontWeight', 'bold',...
        'FontSize', 14)
    ylabel(tt, 'df_n / dlog(' + x_lb + ') (nm^-^1)', 'FontName', 'SansSerif',...
        'FontWeight', 'bold', 'FontSize', 14)
    
    linkaxes([tt1 tt2], 'xy')
    tt1.XLim = 1e9 * [del_d(1), del_d(2)];
    ymax = max([y11; y12; y2]);
    tt1.YLim = [0, ceil(ymax / 0.1) * 0.1];
    
    lgd = legend([p11, p12, p21], {'Preliminary distribution',...
        'Lognormal target', 'Filtered distribution'},...
        'FontName', 'SansSerif', 'FontSize', 12, 'Orientation','horizontal');
    lgd.Layout.Tile = 'north';
end

end

