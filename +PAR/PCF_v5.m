function [g, r] = PCF_v5(pp, nr_o, nr_i, n_g, c_i, c_o, opts)
% "PCF" calculates the Pair Correlation Function (PCF) of a fractal...
%   ...aggregate based on descretizing primary particles and counting...
%   ...the number of grid points in randomly origined radially swept...
%   ...volume elements. 
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pp: Primary particle information matrix
%   n_r1: Number of radial increments determining the resolution of PCF...
%       ...outside the central primary particle.
%   n_r2: Number of radial increments inside the central primary particle.
%   c_o: A coefficient determining the maximum extent of computational...
%       ...domain (this is multiplied by the max pairwise primary...
%       ...particle distance).
%   c_i: A coefficient determining the starting point of radial...
%       ...calculation vector (multipled by min primary particle size).
%   n_g: Number of grid points within the smallest primary particle
%   opts: Options for PCF calculations
% ----------------------------------------------------------------------- %
% 
% Outputs:
%   g: Radial vector of averaged (over primary particles) PCF within and...
%       ...nearby an aggregate.
%   r: The raidal logarithmic points of calculation for PCF starting...
%       ...from inside a central primary particle.
% ----------------------------------------------------------------------- %

% make the options variable if not inputted
if ~exist('opts', 'var') 
    opts = struct();
end

% initialize the visibility variable
if (~isfield(opts, 'vis')) || isempty(opts.vis)
    opts.vis = 'on'; % default to plot the results
end

% initialize the textbar display variable
if (~isfield(opts, 'tbar')) || isempty(opts.tbar)
    opts.tbar = 'on'; % default to print the calculation progress
end

% initialize the logging origin variable
if (~isfield(opts, 'orig')) || isempty(opts.orig)
    opts.orig = 'rand'; % default to start from the center of primaries
end

% initialize figure 
if strcmp(opts.vis, 'on') || strcmp(opts.vis, 'ON') || strcmp(opts.vis, 'On')
    figure;
    h = gcf;
    h.Position = [0, 0, 700, 700];
    set(h, 'color', 'white');
end

% initialize loggingresolution parameters
if ~exist('nr_i', 'var') || isempty(nr_i); nr_i = 20; end
if ~exist('nr_o', 'var') || isempty(nr_o); nr_o = 200; end

% initialize extension parameter
if ~exist('c_i', 'var') || isempty(c_i); c_i = 0.1; end

% initialize extension parameter
if ~exist('c_o', 'var') || isempty(c_o); c_o = 2; end

% initialize grid resolution
if ~exist('n_g', 'var') || isempty(n_g); n_g = 1000; end

n_pp = size(pp,1);
inds_pp = nchoosek(1:n_pp, 2);

% maximum pair-wise primary particle distance
r_o = c_o * max(sqrt(sum((pp(inds_pp(:,1),3:5) - pp(inds_pp(:,2),3:5)).^2, 2)) +...
    (pp(inds_pp(:,1),2) + pp(inds_pp(:,2),2)) / 2);

% the starting point of calculation within the primaries
r_i = c_i * min(pp(:,2));

% the interface bewteen descretization inside and outside the primaries
r_m = geomean(pp(:,2)) / 2;

% initialize pair correlation function vector (along with the vector for...
    % ...radial walk)
g = zeros(nr_i + nr_o + 1, 1);
r = zeros(nr_i + nr_o + 1, 1);

rr_o = (r_o / r_m)^(1 / nr_o); % radial increment factor
r0_o = r_m * ones(nr_o + 1, 1); % initialize auxiliary vector for radial incrementation

% same as above but with non-log increments and for inside the primaries
r0_i = (r_i : (r_m - r_i) / nr_i : r_m)';

% merge the radii inside and outside the primaries
r0 = [0; r0_i(1:end-1); r0_o];

g0 = zeros(1,n_pp); % placeholder for temporary PCF values from each...
    % ...individual primary particle

% make a grid set within the primaries
r_ppdis = cell(n_pp, 1); % placeholder to store grid points
c_ppdis = zeros(n_pp, 1); % initialize lattice concentration within each primary particle
ind_dmin = find(pp(:,2) == min(pp(:,2)),1);
[r_ppdis{ind_dmin}, c_ppdis(ind_dmin)] = PAR.MCDISCRETIZEPP_v2(pp(ind_dmin,2),...
    pp(ind_dmin,3:5), n_g); % discretize smallest primary particle

% discretize the rest of primaries based on the concentration of smallest one
for i = 1 : n_pp
    if i~= ind_dmin
        n_ppdis = round(c_ppdis(ind_dmin) * pi * pp(i,2)^3 / 6);
        [r_ppdis{i}, c_ppdis(i)] = PAR.MCDISCRETIZEPP_v2(pp(i,2), pp(i,3:5), n_ppdis);
    end
end
r_ppdis = cat(1, r_ppdis{:}); % merge the grids made in primaries

% Initialize textbar
if strcmp(opts.tbar, 'on') || strcmp(opts.tbar, 'ON') || strcmp(opts.tbar, 'On')
    fprintf('Volume sweeping started...')
    disp(' ')
    UTILS.TEXTBAR([0, nr_i + nr_o + 1]);
end

for i = 1 : nr_i + nr_o + 1
    % make radial incrementation point
    if i > nr_i + 1
        r0(i+1) = r0(i+1) * rr_o^(i - (nr_i + 1));
        r(i) = sqrt(r0(i+1) * r0(i));
    else
        r(i) = (r0(i+1) + r0(i)) / 2;
    end   
    
    % PCF equals number of grid points within the swept volume over that volume
    for j = 1 : n_pp
        % define the origin of radial volume sweeping
        if strcmp(opts.orig, 'cntr') || strcmp(opts.orig, 'CNTR') || strcmp(opts.orig, 'Cntr')
            rc = pp(j,3:5); % origin to be the center of each primary particle
        elseif strcmp(opts.orig, 'rand') || strcmp(opts.orig, 'RAND') || strcmp(opts.orig, 'Rand')
            % origin to be a random point within each primary particle
            rc0 = rand(1,3); % rho, theta and phi values, respectively, in a cylindrical coordinate
            rc = zeros(1,3); % initialize the origin coordinates
            rc(1) = (pp(j,2) / 2) * rc0(1) .* cos(2 * pi * rc0(2)) .* sin(pi * rc0(3)) + pp(j,3); % x = rho * cos(theta) * sin(phi) 
            rc(2) = (pp(j,2) / 2) * rc0(1) .* sin(2 * pi * rc0(2)) .* sin(pi * rc0(3)) + pp(j,4); % y = rho * sin(theta) * sin(phi) 
            rc(3) = (pp(j,2) / 2) * rc0(1) .* cos(pi * rc0(3)) + pp(j,5); % z = rho * cos(phi)
        end
        
        r_pair = sqrt(sum((r_ppdis - rc).^2, 2));
        
        if i > 1
            g0(j) = nnz((r_pair >= r0(i)) & (r_pair < r0(i+1))) /...
                ((4/3) * pi * (r0(i+1)^3 - r0(i)^3)) / c_ppdis(j);
        else
            g0(j) = nnz(r_pair < r0(i+1)) /...
                ((4/3) * pi * (r0(i+1)^3 - r0(i)^3)) / c_ppdis(j);
        end            
    end
    
    g(i) = mean(g0);
    
    if strcmp(opts.tbar, 'on') || strcmp(opts.tbar, 'ON') || strcmp(opts.tbar, 'On')
        UTILS.TEXTBAR([i, nr_i + nr_o + 1]); % Update textbar
    end
end
% remove noises (partly)
r = r(g~=0);
g = g(g~=0);

r = r / r_m; % normalize the radial increments

% plot PCF vs. normalized radial distance averaged over different primary particles
if strcmp(opts.vis, 'on') || strcmp(opts.vis, 'ON') || strcmp(opts.vis, 'On')
    plot(r(2:end), g(2:end));
    hold on

    box on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18,...
        'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')    
    xlabel('$\overline{r}$ (-)', 'interpreter', 'latex', 'FontSize', 20)
    ylabel('$\overline{g}$($\overline{r}$) (-)', 'interpreter', 'latex', 'FontSize', 20)
end

end

