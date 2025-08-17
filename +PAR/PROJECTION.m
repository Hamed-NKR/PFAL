function [pa_avg, af_avg, h] = PROJECTION(pars, kk_pars, rsl_samp,...
    rsl_avg, kk_h, opts)
% "PROJECTION" computes the orientation averaged projected area of the...
%   ...agggregates based on Monte Carlo (MC) approach.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pars: Particle information structure/class
%   kk_pars: Indices of the aggregates to be analyzed
%   rsl_samp: Sampling resolution (number of spatial points used by MC)...
%       ...at each view angle for each aggregate
%   rsl_avg: Number of random orientations used for averaging the area...
%       ...for each aggregate in MC framework
%   kk_h: ~ to be plotted
% ----------------------------------------------------------------------- %
%
% Outputs:
%   pa_avg: Average projected area array for each aggregate over diffeerent
%       ...angles
%   af_avg: Fraction of orientation averaged projected area to the...
%       ...surrounding rectangle
%   h: Output figure handle (binary response plots of the samples points...
%       ...being inside(1)/outside(0) the aggregate layed over the...
%       ...original aggregate outline)
% ----------------------------------------------------------------------- %

% initialize visibility options variable
if ~exist('opts', 'var') 
    opts = struct();
end
if (~isfield(opts, 'tbar')) || isempty(opts.tbar)
    opts.tbar = 'on'; % textbar on by default
end

% Total number of aggregates
if isa(pars, 'AGG')
    n_tot = length(pars);
else
    n_tot = length(pars.n);
end

% Assigning the aggregates to be analyzed if not given
if ~exist('kk_pars', 'var') || isempty(kk_pars)
    kk_pars = 1 : n_tot;
end

% Determining the sampling resolution if missing
if ~exist('rsl_samp', 'var') || isempty(rsl_samp)
    rsl_samp = 1e4;
elseif rsl_samp < 1e2
    warning('Spatial resolution too low! (better be >= 100)')
end

% Determining the averaging resolution if missing
if ~exist('rsl_avg', 'var') || isempty(rsl_avg)
    rsl_avg = 20;
elseif rsl_avg < 5
    warning('Angular resolution too low! (better be >= 10)')
end
jj = round(1 + (rsl_avg - 1) .* (0 : 1/3 : 1)'); % Indices of angles to...
% ...be plotted

% Initializing the aggregate indices for plotting if not defined
if ~exist('kk_h', 'var')
    kk_h = [];
end

% Initializing the figure handle
if isempty(kk_h)
    h = [];
elseif length(kk_h) > 6
    error('Out of range number of aggregates requested for plotting! (should be <= 6)')
else
    figure;
    h = gcf;
    
    % Setting figure size, position and background
    h.Position = [0, 0, 892.1 * length(kk_h) / 4, 892.1];
    set(h, 'color', 'white');
    
    % Setting the subplots' layout
    tt = tiledlayout(4, length(kk_h));
    tt.TileSpacing = 'compact';
    tt.Padding = 'compact';
    tt.TileIndexing = 'columnmajor';
end

% Compiling the aggregates' properties
if isa(pars, 'AGG')
    pp0 = cat(1, pars(kk_pars).pp);
    n = cat(1, pars(kk_pars).n);
else
    pp0 = cat(1, pars.pp(kk_pars));
    n = cat(1, pars.n(kk_pars));
end

n_agg = length(kk_pars); % Number of aggregates to be analyzed
pa = zeros(n_agg, length(rsl_avg), 2); % Initializing the projected area...
    % ...set (average areas & area fractions)

if strcmp(opts.tbar, 'on') || strcmp(opts.tbar, 'ON') || strcmp(opts.tbar, 'On')
    disp('Computing projected area...');
    UTILS.TEXTBAR([0, n_agg * rsl_avg]); % Initializing 1st layer textbar (for aggergates)
end

for i = 1 : n_agg
%     fprintf('Aggregate %d:', i);
%     disp(newline)
    
    angs = 2 * pi * rand(rsl_avg, 3); % A uniform random set of three...
        % ...intrinsic Euler angles for rotating the aggregtes to...
        % ...average the projected area
    
%     UTILS.TEXTBAR([0, rsl_avg]); % Initializing 2nd layer textbar...
        % ...(for orientations)
    
    for j = 1 : rsl_avg
        pp = cell2mat(PAR.ROTATE(pp0(i), n(i), angs(j,:))); % Rotating...
            % ...the aggregate in random direction and converting to matrix
        
        % Extension range in x & y directions
        rpmax = max(pp(:,2)) / 2;
        x_rng = [min(pp(:,3)) - rpmax, max(pp(:,3)) + rpmax];
        y_rng = [min(pp(:,4)) - rpmax, max(pp(:,4)) + rpmax];
        
        pp_temp = repmat(pp, rsl_samp, 1);
        
        mcpoints = [x_rng(1), y_rng(1)] + [x_rng(2) - x_rng(1),...
            y_rng(2) - y_rng(1)] .* rand(rsl_samp, 2); % Locating the...
                % ...MC uniformly random points 
        
        mcstat = sqrt(sum((pp_temp(:, 3:4) -...
            repelem(mcpoints, n(i), 1)).^2, 2)) <= pp_temp(:, 2) / 2;
                % Check if the MC points fall within primaries
        mcstat = mat2cell(mcstat, n(i) * ones(rsl_samp,1));
        
        % Determine if the MC points are inside the aggregates
        for k = 1 : rsl_samp
            mcstat{k} = any(mcstat{k});
        end
        mcstat = cell2mat(mcstat);
        
        pa(i,j,2) = nnz(mcstat) / length(mcpoints); % Fraction of points...
            % ...falling inside the domain
        
        pa(i,j,1) = abs(pa(i,j,2) * (x_rng(2) - x_rng(1)) *...
            (y_rng(2) - y_rng(1))); % Multiplying by the total domain...
                % ...area to get the projected area
        
        if strcmp(opts.tbar, 'on') || strcmp(opts.tbar, 'ON') || strcmp(opts.tbar, 'On')
            UTILS.TEXTBAR([(i - 1) * rsl_avg + j, n_agg * rsl_avg]) % Updating textbar
        end
        
        % Plotting the results
        if ismember(kk_pars(i), kk_h) && ismember(j, jj) % Plotting only...
                % ...for 4 angles
            nexttile

            % The original outline of aggregate (i)
            pl1 = viscircles([pp(:,3), pp(:,4)], pp(:,2) / 2,...
                'EnhanceVisibility', false, 'Color', [0 0.4470 0.7410],...
                'LineWidth', 1);
            pl1.Children.Color(4) = 0.6;
            hold on
                        
            % MC points
            pl2 = scatter(mcpoints(mcstat,1), mcpoints(mcstat,2), 1,...
                [0.8500 0.3250 0.0980], 'filled');
            
            % Plot specs
            axis equal
            box on
            
            xlim(x_rng)
            ylim(y_rng)
            set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
                'TickLength', [0.04 0.04])
            
            if j == min(jj)
                title({"Agg id: " + num2str(kk_pars(i), "%d"), ""},...
                    'FontName', 'SansSerif', 'FontWeight', 'bold',...
                    'FontSize', 14);
            end
            
            if find(kk_h == kk_pars(i)) == 1 &&...
                    find(jj == j) == 1 % Only one legend for the whole plot
                lgd = legend([pl1, pl2],...
                {'Agg outline', 'MC points'},'Orientation','horizontal',...
                'FontName', 'SansSerif', 'FontSize', 12);
                lgd.Layout.Tile = 'north';
            end
            
%             UTILS.TEXTBAR([j, rsl_avg]); % Completing orientation j
        end
    end
%     disp(newline)
%     UTILS.TEXTBAR([i, n_agg]); % Finishing j
%     disp(newline)
end

if ~isempty(kk_h)
    title(tt, 'Aggregate projected area in different orientations',...
        'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 18)
    subtitle(tt, ' ','FontName', 'SansSerif', 'FontSize', 8)
    xlabel(tt, {'\fontsize{4} ', '\fontsize{14}x (m)'},'interpreter','tex',...
        'FontName', 'SansSerif', 'FontWeight', 'bold')
    ylabel(tt, {'\fontsize{14}y (m)', '\fontsize{4} '},'interpreter','tex',...
        'FontName', 'SansSerif', 'FontWeight', 'bold')
end

pa_avg = mean(pa(:,:,1), 2); % Averaging are over orientations
af_avg = mean(pa(:,:,2), 2); % Averaged area fraction

if nargout < 3
    clear h; % Deleting figure handle if not requested as an output
end

end

