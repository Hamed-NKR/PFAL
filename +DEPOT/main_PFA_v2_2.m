clc
clear
clf('reset')
close all

%% Initializations %%

[params_ud, params_const] = TRANSP.INIT_PARAMS('MCEM_PFAParams'); % Read the input file

% Stage 1:
n_stor = 5; % Number of data storage occurrences
n_try = 3; % Number of DLCA trials

gstd_dppi_ens = 1.4; % Geometric standard deviation of ensemble average primary particle size

% The desired distribution parameters for aggregates after lognormal sampling
mu_da_glob = 1.25e-7;
std_da_glob = 1.4;

npp_min = 10; % Aggregate filtering criterion
npp_max = 200; % Iteration limit parameter in terms of number of primaries within the aggregate

j_max = 1e7; % Stage 1 marching index limit

opts.visual = 'on'; % flage for display of lognormal sampling process
opts.randvar = 'area'; % flag for type of size used in lognormal sampling

% flag for caclculation method of mobility diameter
opts_mobil.mtd = 'interp';
% opts_mobil.mtd = 'continuum';

% assign marching timecale multiplier
% opts_mobil.c_dt = 0.1;

% flag for collision type determination if applicable
% opts_grow.col = 'coal';

% Stage 2
f_dil = 0.1; % Dilution factor for post-flame agglomeration

D_TEM = 0.35; % polydispersity exponent

c_proj = 3; % the coefficient to improve the resolution of projected area...
    % ...calculation in the second stage

k_max = 1e5; % Iteration limit parameter

n_kk = 5; % Number of saving timespots (in round2)

opts2.plotlabel = 'on'; % flag to display lablels of gstd on dp vs da plots
opts2.ploteb = 'on'; % error bar display flag

opts2.savelast = 'on'; % flag to whether save the last piece of data

opts2.datastore = 'n_agg'; % flag to decide on criterion for data saving...
    % ...('n_hyb' for number of hybridity regions, 'n_agg' for total...
    % ...number of aggs within the domain)

% Growth limit extremes
if strcmp(opts2.datastore, 'n_hyb')
    kk_min = 2;
    kk_max = 10;
else
    kk_min = 0.5;
    kk_max = 0.02;
end

opts2_kin.visual = 'on'; % flag to visualization of kinetic properties

%% 1st stage %%

% params_ud.Value(1) = 1e-3; % reassign volume fraction (don't use unless for multiple serial runs)

pp0 = cell(n_stor, n_try); % Primary particle data storage cell array for initial monodisperse aggregation

pp0_n = cell(n_stor, n_try); % Placeholder for number of primaries withing aggs

i_pp0 = cell(n_stor,1); % Placeholder for index of primary particles

% Iteration steps to be stored
r_str = (npp_max / npp_min)^(1 / (n_stor - 1));
jj = zeros(n_stor,1);
for i = 1 : n_stor
    jj(i) = round(npp_min * (r_str^(i-1)));
end

fprintf('Simulating monodisperse clustering...')
disp(newline)

[pars, fl] = TRANSP.INIT_DOM(params_ud, params_const); % Initialize particle and fluid structs

% Make an initial lognormal distibution of ensemble average size of primaries for classic dlca trials
dppi = lognrnd(log(params_ud.Value(8)), log(gstd_dppi_ens), [n_try,1]);

% Make a placeholder for diffusive properties
mbl_stor = cell(n_stor, n_try); % the particle-resolved values for selected moments 
mbl_rt.n = zeros(j_max, n_try); % real-time total number of aggregates
mbl_rt.t = zeros(j_max, n_try); % time at each iteration
mbl_rt.delt = zeros(j_max, 2 * n_try); % real-time mean and std. of calculated timesteps of aggregates 
mbl_rt.tau = zeros(j_max, 2 * n_try); % real-time mean and std. of characteristic time-scales of aggregates
mbl_rt.kn_kin = zeros(j_max, 2 * n_try); % real-time mean and std. of kinetic Knudsen number
mbl_rt.kn_diff = zeros(j_max, 2 * n_try); % real-time mean and std. of diffusive Knudsen number

% Generate a library of monodisperse aggregates with classic DLCA
for i = 1 : n_try
    fprintf('trial %d:', i)
    disp(' ')
    
    % Initilize monodisperse pps for classic DLCA
    [pp_d, pars.n] = PAR.INIT_DIAM(params_ud.Value(5), params_ud.Value(6:7),...
        [dppi(i); params_ud.Value(9:10)]); % Initialize pp sizes

    pars.pp = mat2cell([(1:size(pp_d))', pp_d, zeros(size(pp_d,1),3),...
        (1:size(pp_d))'], pars.n); % Assign pp indices and sizes

    if params_ud.Value(6) ~= 0
        pars.pp = PAR.INIT_MORPH_RAND(pars.pp); % Randomly initialize pp locations within aggs
    end

    [pars, params_ud] = PAR.INIT_LOC(pars, params_ud); % Randomize initial locations

    pars = PAR.SIZING(pars); % Calculate sizes

    pars = TRANSP.MOBIL(pars, fl, params_const, opts_mobil); % Calculate mobility props

    % record real-time kinetics-related data
    mbl_rt.n(1,i) = length(pars.n);
    mbl_rt.t(1,i) = min(pars.delt);
    mbl_rt.delt(1,2*(i-1)+(1:2)) = [mean(pars.delt), std(pars.delt)];
    mbl_rt.tau(1,2*(i-1)+(1:2)) = [mean(pars.tau), std(pars.tau)];
    mbl_rt.kn_kin(1,2*(i-1)+(1:2)) = [mean(pars.kn_kin), std(pars.kn_kin)];
    mbl_rt.kn_diff(1,2*(i-1)+(1:2)) = [mean(pars.kn_diff), std(pars.kn_diff)];

    pars.v = PAR.INIT_VEL(pars.pp, pars.n, fl.temp, params_const); % Randomize initial velocities

    j = 2; % initialize iteration index
    jjj = 1; % pp data storage tracking index
    
    fprintf('\b')
    fprintf('DLCA simulation...')
    disp(' ')
    UTILS.TEXTBAR([0, j_max]); % Initialize textbar
    UTILS.TEXTBAR([1, j_max]); % Iteration 1 already done
    
    % 1st stage DLCA iterations
    while (j <= j_max) &&...
            (sum(cat(1, pars.n) >= npp_max) < round(0.9 * length(cat(1, pars.n)))) &&...
            (length(pars.n) > 1)
        pars = TRANSP.MARCH(pars, fl, params_const); % Solve for transport
        
        pars = TRANSP.PBC(params_ud.Value(2:4), pars); % Apply periodic BCs
        
        pars = COL.GROW(pars); % Cluster the particles
        
        if sum(cat(1, pars.n) >= jj(jjj)) >= round(0.9 * length(cat(1, pars.n)))
            
            % Store number of primaries
            pp0_n{jjj,i} = cat(1, pars.n);
            
            %         % Remove small aggregates
            %         if jjj == 1
            %             pars.pp(pp0_n{jjj} < npp_min) = [];
            %             pars.n(pp0_n{jjj} < npp_min) = [];
            %             pars.r(pp0_n{jjj} < npp_min, :) = [];
            %             pars.v(pp0_n{jjj} < npp_min, :) = [];
            %             pp0_n{jjj}(pp0_n{jjj} < npp_min) = [];
            %         end
            
            pp0{jjj,i} = pars.pp; % Store pp info
            mbl_stor{jjj,i}.kn_diff = pars.kn_diff; % diffusive Knudsen no.
            mbl_stor{jjj,i}.kn_kin = pars.kn_kin; % kinetic Knudsen no.
            
            % Update pp indices
            i_pp0{jjj} = cell(length(cat(1, pars.n)), 1);
            for j4 = 1 : length(cat(1, pars.n))
                i_pp0{jjj}{j4} = pars.pp{j4}(:,1);
                if (i > 1) || (jjj > 1)
                    pp0{jjj,i}{j4}(:,1) = pp0{jjj,i}{j4}(:,1) +...
                        ((i - 1) * n_stor + jjj - 1) * params_ud.Value(5);
                end
            end
            
            jjj = jjj + 1;
        end
                
        pars = PAR.SIZING(pars); % Update sizes
        
        pars = TRANSP.MOBIL(pars, fl, params_const, opts_mobil); % Update mobility properties
        
        % record real-time kinetics-related data
        mbl_rt.n(j+1,i) = length(pars.n);
        mbl_rt.t(j+1,i) = min(pars.delt) + mbl_rt.t(j,i);
        mbl_rt.delt(j+1,2*(i-1)+(1:2)) = [mean(pars.delt), std(pars.delt)];
        mbl_rt.tau(j+1,2*(i-1)+(1:2)) = [mean(pars.tau), std(pars.tau)];
        mbl_rt.kn_kin(j+1,2*(i-1)+(1:2)) = [mean(pars.kn_kin), std(pars.kn_kin)];
        mbl_rt.kn_diff(j+1,2*(i-1)+(1:2)) = [mean(pars.kn_diff), std(pars.kn_diff)];
        
        UTILS.TEXTBAR([j, j_max]); % Update textbar
        
        j = j + 1; % Update iteration ind.
    end
    
    % Store the pp info and number of pps and update the indices for the last step
    if jjj <= n_stor
        pp0{jjj,i} = pars.pp;
        pp0_n{jjj,i} = cat(1, pars.n);
        i_pp0{jjj} = cell(length(cat(1, pars.n)), 1);
        for j4 = 1 : length(cat(1, pars.n))
            i_pp0{jjj}{j4} = pars.pp{j4}(:,1);
            if (i > 1) || (jjj > 1)
                pp0{jjj,i}{j4}(:,1) = pp0{jjj,i}{j4}(:,1) +...
                    ((i - 1) * n_stor + jjj - 1) * params_ud.Value(5);
            end
        end
    end
    
    % Remove similar aggs
    if jjj > 1
        i_pp00 = cat(1, i_pp0{:}); % Compile indices
        ij = nchoosek(1 : length(cell2mat(pp0_n(:,i))), 2);
        ind_rmv = [];
        for ii = 1 : length(ij)
            if isequal(mod(i_pp00{ij(ii,1)}(:,1), params_ud.Value(5)),...
                    mod(i_pp00{ij(ii,2)}(:,1), params_ud.Value(5)))
                ind_rmv = [ind_rmv; ij(ii,2)];
            end
        end
        
        if ~isempty(ind_rmv)
            ind_rmv = unique(ind_rmv);
            n_rmv = length(ind_rmv);
            
            agg0_n = cell2mat(cellfun(@size, pp0(:,i), 'UniformOutput', false));
            agg0_n(:,2) = [];
            agg0_n_c = [0; cumsum(agg0_n)];
            
            ind_rmv1 = zeros(n_rmv,1);
            ind_rmv2 = zeros(n_rmv,1);
            for iii = 1 : n_rmv
                ind_rmv1(iii) = find(ind_rmv(iii) <= agg0_n_c, 1) - 1;
                ind_rmv2(iii) = ind_rmv(iii) - agg0_n_c(ind_rmv1(iii));
            end
            
            ind_rmv11 = unique(ind_rmv1);
            n_rmv1 = length(ind_rmv11);
            for iii = 1 : n_rmv1
                i4 = find(ind_rmv1 == ind_rmv11(iii));
                pp0{ind_rmv11(iii),i}(ind_rmv2(i4)) = [];
                pp0_n{ind_rmv11(iii),i}(ind_rmv2(i4)) = [];
            end
        end
    end
    
    disp(newline)
end

% Convert the pp data into a 1d cell array
pp0 = pp0(:);
pp0_n = pp0_n(:);

% Remove unused cells
pp0 = pp0(~cellfun('isempty', pp0));
pp0_n = pp0_n(~cellfun('isempty', pp0_n));

pp0 = cat(1, pp0{:}); % Merge pp info from different times
pp0_n = cat(1, pp0_n{:}); % Compile number of primaries data

n_pp0 = sum(pp0_n); % total number of primary particles after stage 1
n_agg0 = size(pp0, 1); % % total number of aggregates before rescaling

dpp0_g = PAR.GEOMEANPP(pp0);
dpp0_g = dpp0_g(:,1); % Current mean pp size withing aggs
pp00 = pp0;
dpp_emh = ((17.8^(1 / D_TEM) / 100) * (pp0_n / 1.1).^...
    (1 / (2 * 1.08))).^(D_TEM / (1 - D_TEM)); % Desired mean pp size based on external mixing hypothesis
r_dpp = 1e-9 * dpp_emh ./ dpp0_g; % Size conversion ratio
% Rescale stored aggregates
for i = 1 : n_agg0
    pp00{i}(:,2:5) = pp00{i}(:,2:5) * r_dpp(i); % Rescale primary particle size
end

pars = structfun(@(x) [], pars, 'UniformOutput', false); % Clear pars struct fields
pars.pp = pp00; % Assign pp info for projected area calculations
pars.n = pp0_n;

da0 = 2 * sqrt(PAR.PROJECTION(pars, [], 1e4, 20) / pi); % Get projected...
    % ...area diameter for monodisperse populations

dpp00_g = PAR.GEOMEANPP(pars.pp);
dpp00_g = dpp00_g(:,1); % Mean primary particle diameter

% Filter particles for a lognormal target area distribution
opts.nfit = 'on';
[da1, ind1] = TRANSP.LNSAMPLING(da0, mu_da_glob, std_da_glob, [], 20,...
    [1e-8, 1e-6], opts);

da1 = da1(~cellfun('isempty', da1)); % Remove empty bins
da1 = cat(1, da1{:}); % Compile the data

ind1 = ind1(~cellfun('isempty', ind1));
ind1 = cat(1, ind1{:});

dpp1_g = dpp00_g(ind1);
pp1 = pp00(ind1);
pp1_n = pp0_n(ind1);

n_agg1 = size(pp1, 1);
pp11 = cat(1, pp1{:});
gstd_dpp1_ens = UTILS.GEOSTD(pp11(:,2)); % ensemble average pp standard deviation after lognormal sampling

% Update the indices of particles due to existence of duplicates
i_pp1 = zeros(n_agg1,1);
for i = 1 : n_agg1
    i_pp1(i) = max(pp1{i}(:,1));
end
[~, i_pp11] = sort(i_pp1);
dpp1_g = dpp1_g(i_pp11);
pp1 = pp1(i_pp11);
pp1_n = pp1_n(i_pp11);
da1 = da1(i_pp11);
pp1_temp = cell2mat(pp1);
pp1_temp(:,1) = pp1_temp(:,1) +...
    repelem(params_ud.Value(5) * (0 : n_agg1 - 1)', pp1_n) +...
    n_stor * n_try * params_ud.Value(5);
pp1 = mat2cell(pp1_temp, pp1_n);

% Randomize order of aggregates
i_rnd = randperm(n_agg1);

dpp1_g = dpp1_g(i_rnd);
pp1 = pp1(i_rnd);
pp1_n = pp1_n(i_rnd);
da1 = da1(i_rnd);

pars.pp = pp1; % Assign pp info for stage 2
pars.n = pp1_n; % Assign number distribution of primaries

%% 2nd stage %%

d1max = PAR.TERRITORY(pp1, pp1_n);
pp11 = cat(1, pp1{:});
r_vfadj = sum(d1max .^ 3) / sum(pp11(:,2) .^ 3); % Adjust for real vs. max. extent volume fraction
params_ud.Value(1) = params_ud.Value(1) * f_dil * r_vfadj; % Dilute the concentration

[pars, params_ud] = PAR.INIT_LOC(pars, params_ud); % Assign random locations to aggregates

pars = PAR.SIZING(pars); % Get sizes

% change the gas properties to room condition
opts_fl.amb = 'room';
[fl.mu, fl.lambda] = TRANSP.FLPROPS(fl, params_const, opts_fl);

pars = TRANSP.MOBIL(pars, fl, params_const, opts_mobil); % Get mobility props

pars.v = PAR.INIT_VEL(pars.pp, pars.n, fl.temp, params_const); % Assign random velocities to aggregates

opts2.indupdate = 'off'; % flag to update aggregate ids upon each collision
    % ...(disabled to track hybridty)

k = 2; % Initialize iteration index

time = zeros(k_max,1); % Initialize time storage array
n_agg2 = zeros(k_max,1); % Real time array of aggregate counts during post-flame agglomeration
n_agg2(1) = length(pars.n);

if n_kk > 10
    n_kk = 10; % Set a limit for data storage
end

% Generate data saving timespots
% kk = kk_min : (kk_max - kk_min) / (n_kk - 1) : kk_max;
r_kk = (kk_max / kk_min)^(1 / (n_kk - 1));
kk = kk_min * ones(n_kk,1);
for i = 2 : n_kk
    kk(i) = kk(i) * r_kk^(i-1);
end
if (strcmp(opts2.datastore, 'n_hyb'))
    kk = unique(round(kk)); % Round the values generated for number of regions
    n_kk = length(kk);
end

kkk = 1; % Index for data saving timespots

pars.n_hyb = ones(n_agg1,1); % All aggs are initially near-monodisperse

parsdata = struct('dpp', cell(n_kk + 1, 1), 'dpp_g', cell(n_kk + 1, 1),...
    'da', cell(n_kk + 1, 1), 'dg', cell(n_kk + 1, 1), 'npp', cell(n_kk + 1, 1),...
    'n_hyb', cell(n_kk + 1, 1), 'pp', cell(n_kk + 1, 1)); % Placeholder for hybrid aggregates data

% Assign values to the dataset from the already available monodipserse population
parsdata(1).da = da1;
parsdata(1).dpp = PAR.MEANPP(pp1);
parsdata(1).dpp_g = PAR.GEOMEANPP(pp1);
parsdata(1).dg = PAR.GYRATION(pp1, pp1_n);
parsdata(1).npp = pp1_n;

disp(' ')
disp('Simulating post-flame agglomeration...')
UTILS.TEXTBAR([0, k_max]); % Initialize textbar
UTILS.TEXTBAR([1, k_max]); % Iteration 1 already done

% Produce hybridized particles with post-flame agglomeration
while (k <= k_max) && (kkk <= n_kk) && (length(pars.n) > 1)
    if (strcmp(opts2.datastore, 'n_hyb') &&...
            (nnz(pars.n_hyb >= kk(kkk)) / n_agg2(k-1) >= 0.9)) ||...
            (strcmp(opts2.datastore, 'n_agg') &&...
            (n_agg2(k-1) <= kk(kkk) * n_agg2(1)))
        
        % save aggregate data over time
        parsdata(kkk + 1).dpp = pars.dpp;
        parsdata(kkk + 1).dpp_g = pars.dpp_g;
        parsdata(kkk + 1).dg = pars.dg;
        parsdata(kkk + 1).npp = pars.n;
        parsdata(kkk + 1).n_hyb = pars.n_hyb;
        parsdata(kkk + 1).pp = pars.pp;
        
        disp(' ')
        parsdata(kkk + 1).da = 2 * sqrt(PAR.PROJECTION(pars, [],...
            min(3e4, max(c_proj * 1e4, 10 * max(pars.n))), 20) / pi); % Get projected area size
        disp(' ')
                
        kkk = kkk + 1; % Go to next saving spot
    end
    
    [pars, delt] = TRANSP.MARCH(pars, fl, params_const); % Solve for transport
    
    pars = TRANSP.PBC(params_ud.Value(2:4), pars); % Apply periodic BCs
    
    pars = COL.GROW(pars, opts2); % Cluster the particles
    
    pars.n_hyb = COL.HYBRIDITY(pars.pp, pars.n); % Count the number of...
        % ...monodisperse regions within a hybrid polydisperse...
        % ...aggregates formed by post-flame agglomeration
    
    pars = PAR.SIZING(pars); % Update sizes
    
    pars = TRANSP.MOBIL(pars, fl, params_const, opts_mobil); % Update mobility properties
    
    time(k) = time(k-1) + delt; % Update time
    n_agg2(k) = length(pars.n); % Record number of aggs
    
    UTILS.TEXTBAR([k, k_max]); % Update textbar
    
    k = k + 1; % Update iteration ind.
end

% Save last set of pp size data if requested
if (kkk <= n_kk) && (strcmp(opts2.savelast, 'on'))
    parsdata(kkk + 1).dpp = pars.dpp;
    parsdata(kkk + 1).dpp_g = pars.dpp_g;
    parsdata(kkk + 1).dg = pars.dg;
    parsdata(kkk + 1).npp = pars.n;
    parsdata(kkk + 1).da = 2 * sqrt(PAR.PROJECTION(pars, [],...
        min(3e4, max(c_proj * 1e4, 10 * max(pars.n))), 20) / pi);
    
    kkk = kkk + 1;
end

% Remove the unused cells and update the dataset number
if kkk <= n_kk
    parsdata(kkk + 1 : end).dpp = [];
    parsdata(kkk + 1 : end).dpp_g = [];
    parsdata(kkk + 1 : end).dg = [];
    parsdata(kkk + 1 : end).npp = [];
    parsdata(kkk + 1 : end).da = [];
    
    n_kk = kkk - 1;
end

% Initialize dp vs da figure
figure(2)
h2 = gcf;
if ~all(h2.Position == [0, 0, 600, 600])
    h2.Position = [0, 0, 600, 600];
end
set(h2, 'color', 'white')

dpp_uc = linspace(5, 65, 1000);
da_uc = 100 * (dpp_uc / 17.8).^(1 / 0.35);

p2_1 = plot(da_uc, dpp_uc, 'Color', [0.5 0.5 0.5],...
    'LineStyle', '-.', 'LineWidth', 2.5); % Plot universal correlation
hold on

p2_2 = scatter(1e9 * da1, 1e9 * dpp1_g, 25, [0.5 0.5 0.5], 'o'); % Plot monodisperse aggs

p2_3 = cell(n_kk,1); % Initialize plot variable for lifetime of hybrids
lgd2_3 = cell(n_kk,1); % Legend text placeholder

ms = [25, 35, 30, 25, 35, 25, 25, 25, 25, 25]; % Marker sizes

mc = colormap(jet); % Marker colormap for visualization of hybrids over aging
ii = round(1 + (length(mc) - 1) .* (0.05 : 0.9 / (n_kk - 1) : 0.95)'); % Descritize the map based on number of life stages
mc = mc(ii,:); % Get the descretized colormap
mc = flip(mc,1); % Reverse the map direction

mt = {'^', 's', 'd', 'v', 'h', '+', '<', 'p', 'x', '>'}; % Marker type depot

for i = 1 : n_kk
    p2_3{i} = scatter(1e9 * parsdata(i+1).da, 1e9 * parsdata(i+1).dpp_g(:,1),...
        ms(i), mc(i,:), 'filled', mt{i}); % Plot hybrid aggs
    
    if strcmp(opts2.plotlabel, 'on')
        lbl = num2str(parsdata(i+1).dpp_g(:,2), '%.2f');
        text(1e9 * parsdata(i+1).da, 1e9 * parsdata(i+1).dpp_g(:,1),...
            lbl, 'VerticalAlignment', 'top',...
            'HorizontalAlignment', 'left', 'FontSize', 10)
    end
    
    if strcmp(opts2.ploteb, 'on')
        e_p = parsdata(i+1).dpp_g(:,1) .* abs(parsdata(i+1).dpp_g(:,2) - 1);
        e_n = parsdata(i+1).dpp_g(:,1) .* abs(1 - 1 ./ parsdata(i+1).dpp_g(:,2));
        eb = errorbar(1e9 * parsdata(i+1).da, 1e9 * parsdata(i+1).dpp_g(:,1),...
            1e9 * e_n, 1e9 * e_p, '.');
        eb.Color = mc(i,:);
    end
    
    if strcmp(opts2.datastore, 'n_hyb')
        if i == 1
%             lgd2_3{i} = strcat('Hybrid - \tau = ', num2str(i, '%.1e'),...
%                 ' (s) - n^9^0^p = ', num2str(kk(i)));
            lgd2_3{i} = strcat('Hybrid, n_{hyb}^{90p} =', {' '}, num2str(kk(i), '%d'));
        else
%             lgd2_3{i} = strcat('~ - \tau = ', num2str(i, '%.1e'),...
%                 ' (s) - n^9^0^p = ', num2str(kk(i)));
            lgd2_3{i} = strcat('~, n_{hyb}^{90p} =', {' '}, num2str(kk(i), '%d'));
        end
    else
        if i == 1
%             lgd2_3{i} = strcat('$Hybrid, t =', {' '}, num2str(i, '%.1e'),...
%                 ' {\tau}, {\sigma}_{g,d__{pp}}^{90p} =', {' '}, num2str(kk(i)), '$');
            lgd2_3{i} = strcat('Hybrid, n_{agg} / n_{agg_0} =', {' '}, num2str(kk(i), '%.2f'));
        else
%             lgd2_3{i} = strcat('$~, t =', {' '}, num2str(i, '%.1e'),...
%                 ' {\tau}, {\sigma}_{g,d__{pp}}^{90p} =', {' '}, num2str(kk(i)), '$');
            lgd2_3{i} = strcat('~, n_{agg} / n_{agg_0} =', {' '}, num2str(kk(i), '%.2f'));
        end
    end
end

% axis equal
box on
set(gca, 'FontName', 'SansSerif', 'FontSize', 12, 'TickLength', [0.02 0.02])
xlabel('d_a (nm)', 'FontName', 'SansSerif', 'FontWeight', 'bold',...
    'FontSize', 14)
ylabel('d_p_p (nm)', 'FontName', 'SansSerif', 'FontWeight', 'bold',...
    'FontSize', 14)
ylim([5, 65])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
legend([p2_1, p2_2, cat(1, p2_3{:})'], cat(2,{'Universal correlation',...
    'Monodisperse'}, lgd2_3{:}),...
    'Location', 'northwest', 'FontName', 'SansSerif', 'FontSize', 12);
title('Primary particle size vs projected area equivalent size',...
    'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 16)

% Remove empty cells from number/time data
if k - 1 < k_max
    n_agg2(k : end) = [];
    time(k : end) = [];
end
[beta, tau, z] = TRANSP.KINETIC(n_agg2, time, opts2_kin); % Compute the kinetic properties of aggregation

% figure(4)
% UTILS.RENDER(pars); % display final aggregates

