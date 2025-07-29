%%% This script is still buggy, use the first half of...
    % ...script: main_PFA_v2_2.m instead to generate a library of DLCA...
    % ...aggregates.

clc
clear
clf('reset')
close all

%% initialize simulation variables %%

% read the batch file containing user's inputs
[params0_ud, params0_const] = TRANSP.INIT_PARAMS('LD1_Params');

n_temporal = 10; % set number of steps that data is stored during each...
% ...simulation
n_trial = 3; % set number of separate Langevin dynamics (LD) simulation...
% ...trials to generate library

% geometric standard deviation of ensemble geometric mean of...
% ...primary particle size across LD simulations (not a critical...
% ...parameter due to the scaling imposed later)
s_pp0 = 1.4;

npp0_min = 10; % number of primary particles for smallest aggregates to be stored
npp0_max = 1000; % number of primary particles for largest aggregates to be stored

j_max = 3e4; % maximum number of iterations in each LD simulation

% flag for caclculation method of mobility diameter
% opts0_mobil.mtd = 'interp';
opts0_mobil.mtd = 'continuum';

% turn off progress bar for projected area calculation
opts0_prj.tbar = 'off';

% assign a multiplier to the timesteps of random walk
opts0_mobil.c_dt = 1000;

% flag for whether collision results in coalescence or agglomeration
% opts_grow.col = 'coal';

% resolution of Monte Carlo method for projected area calculation
n0_mc_prj = 1e2;
n0_ang_prj = 5;

%% set up the particle properties for the start of simulations %%

% make a 2d cell array to store primary particle data over time and for
% ...different simulations
pp0 = cell(n_temporal, n_trial);

% placeholder for number of primary particle withing aggregates stored...
% ...in pp0
pp0_n = cell(n_temporal, n_trial);

% placeholder for index of primary particles
i_pp0 = cell(n_temporal,1);

% find number of primary particle within aggregates that correspond to...
% ...intermediate iteration data to be saved
r_temporal = (npp0_max / npp0_min)^(1 / (n_temporal - 1));
npp0_temporal = zeros(n_temporal,1);
for i = 1 : n_temporal
    npp0_temporal(i) = round(npp0_min * (r_temporal^(i-1)));
end

fprintf('initializing simulation of first-stage aggregation (to generate data library)...')
disp(newline)

% initialize particle and fluid structures
[pars0, fl0] = TRANSP.INIT_DOM(params0_ud, params0_const);

% generate an lognormal distibution to select the ensemble geometric...
% ...mean of primary particle diameter for each simulation trial
dpp0_ens = lognrnd(log(params0_ud.Value(8)), log(s_pp0), [n_trial,1]);

% allocate the placeholder for saving average values for properties of...
% ...interest at every timestep

ensdata0.n_agg = zeros(j_max, n_trial); % total number of aggregates
ensdata0.t = zeros(j_max, n_trial); % time at each iteration

% mean and std. of characteristic time-scales of aggregates
ensdata0.tau = zeros(j_max, 2 * n_trial);

% mean and std. of kinetic Knudsen number
ensdata0.kn_kin = zeros(j_max, 2 * n_trial);

% mean and std. of diffusive Knudsen number
ensdata0.kn_diff = zeros(j_max, 2 * n_trial);

% geometric mean and std. of primary particle diameter
ensdata0.dpp = zeros(j_max, 2 * n_trial);

% geometric mean and std. of geometric standard deviation of primary particle diameter
ensdata0.sigmapp = zeros(j_max, 2 * n_trial);

% geometric mean and std. of projected area diameter
ensdata0.da = zeros(j_max, 2 * n_trial);

% geometric mean and std. of mobility diameter
ensdata0.dm = zeros(j_max, 2 * n_trial);

% allocate the placeholder to save average values for properties of...
% ...interest at every timestep
parsdata0 = struct('dpp', cell(n_temporal, n_trial),...
    'sigmapp', cell(n_temporal, n_trial), 'da', cell(n_temporal, n_trial),...
    'dm', cell(n_temporal, n_trial), 'dg', cell(n_temporal, n_trial));

%% start the iterations %%

% nested loop for generating and compiling a library of aggregates...
% ...supposedly forming within the flame
for i = 1 : n_trial % run n_trial number of Langevin dynamics (LD) simulations
    fprintf('trial %d:', i)
    disp(' ')

    % generate an initial distribution of particle sizes
    [pp0_d, pars0.n] = PAR.INIT_DIAM(params0_ud.Value(5), params0_ud.Value(6:7),...
        [dpp0_ens(i); params0_ud.Value(9:10)]); % Initialize pp sizes

    % make the particle structure to store real-time data of aggregate...
    % ...population
    pars0.pp = mat2cell([(1:size(pp0_d))', pp0_d, zeros(size(pp0_d,1),3),...
        (1:size(pp0_d))'], pars0.n);

    % if starting with aggregates rather than primary particles,...
    % ...make a population of uniform random aggregates (non-DLCA)
    if params0_ud.Value(6) ~= 0
        pars0.pp = PAR.INIT_MORPH_RAND(pars0.pp);
    end

    % assign random initial locations
    [pars0, params0_ud] = PAR.INIT_LOC(pars0, params0_ud);

    % calculate initial characteristic sizes
    pars0 = PAR.SIZING(pars0);

    % calculate initial mobility properties
    pars0 = TRANSP.MOBIL(pars0, fl0, params0_const, opts0_mobil);

    % record ensemble-averaged properties at the beginning of simulations
    ensdata0.t(1,i) = min(pars0.delt);
    ensdata0.n_agg(1,i) = length(pars0.n);
    ensdata0.tau(1,2*(i-1)+(1:2)) = [mean(pars0.tau), std(pars0.tau)];
    ensdata0.kn_kin(1,2*(i-1)+(1:2)) = [mean(pars0.kn_kin), std(pars0.kn_kin)];
    ensdata0.kn_diff(1,2*(i-1)+(1:2)) = [mean(pars0.kn_diff), std(pars0.kn_diff)];
    ensdata0.dpp(1,2*(i-1)+(1:2)) = [geomean(pars0.dpp_g(:,1)),...
        UTILS.GEOSTD(pars0.dpp_g(:,1))];
    ensdata0.sigmapp(1,2*(i-1)+(1:2)) = [geomean(pars0.dpp_g(:,2)),...
        UTILS.GEOSTD(pars0.dpp_g(:,2))];
    ensdata0.dm(1,2*(i-1)+(1:2)) = [geomean(pars0.dm), UTILS.GEOSTD(pars0.dm)];
    if strcmp(opts0_mobil.mtd, 'interp')
        ensdata0.da(1,2*(i-1)+(1:2)) = [geomean(pars0.da), UTILS.GEOSTD(pars0.da)];
    end

    % assign random (Maxwell Boltzmann) initial velocities
    pars0.v = PAR.INIT_VEL(pars0.pp, pars0.n, fl0.temp, params0_const);

    j = 2; % initialize iteration index

    jj = 1; % index for primary particle particle data storage tracking

    fprintf('\b')
    fprintf('Simulating...')
    disp(' ')
    UTILS.TEXTBAR([0, j_max]); % initialize textbar
    UTILS.TEXTBAR([1, j_max]); % first iteration already completed

    % iterations starting
    while (j <= j_max) && (sum(cat(1, pars0.n) >= npp0_max) <...
            round(0.9 * length(cat(1, pars0.n)))) && (length(pars0.n) > 1)
        % stop when reaching j_max, or when...
        % ...90% of aggregates reach npp0_max, or when...
        % ...only one aggregate remains

        % solve transport equations and move particles
        pars0 = TRANSP.MARCH(pars0, fl0, params0_const);

        % apply periodic boundary conditions
        pars0 = TRANSP.PBC(params0_ud.Value(2:4), pars0);

        % connect colliding particles
        pars0 = COL.GROW(pars0);

        % update characteristic sizes
        pars0 = PAR.SIZING(pars0);
        
        % update mobility properties
        pars0 = TRANSP.MOBIL(pars0, fl0, params0_const, opts0_mobil);        

        if sum(cat(1, pars0.n) >= npp0_temporal(jj)) >=...
                round(0.9 * length(cat(1, pars0.n)))
            % save data when 90% of aggregates reach npp0_temporal(jj)

            % save number of primary particles (this is needed later to...
            % ...remove duplicates)
            pp0_n{jj,i} = cat(1, pars0.n);

            % % Remove aggregates smaller than npp_min
            % if jj == 1
            %     pars.pp(pp0_n{jj} < npp_min) = [];
            %     pars.n(pp0_n{jj} < npp_min) = [];
            %     pars.r(pp0_n{jj} < npp_min, :) = [];
            %     pars.v(pp0_n{jj} < npp_min, :) = [];
            %     pp0_n{jj}(pp0_n{jj} < npp_min) = [];
            % end
            
            % save primary particle sizes and coordinates
            pp0{jj,i} = pars0.pp;
            
            % calculate projected area diameter if not calculated already
            if ~strcmp(opts0_mobil.mtd, 'interp')
                pars0.da = 2 * sqrt(PAR.PROJECTION(pars0, [],...
                    n0_mc_prj, n0_ang_prj, [], opts0_prj) / pi);
            end

            % save properties of aggregates
            parsdata0(jj,i).dpp = pars0.dpp_g(:,1);
            parsdata0(jj,i).sigmapp = pars0.dpp_g(:,2);
            parsdata0(jj,i).dm = pars0.dm;
            parsdata0(jj,i).dg = pars0.dg;
            parsdata0(jj,i).da = pars0.da;

            % update primary particle indices (needed to avoid confusion...
            % ...when identifying duplicates)
            i_pp0{jj} = cell(length(cat(1, pars0.n)), 1);
            for jjj = 1 : length(cat(1, pars0.n))
                i_pp0{jj}{jjj} = pars0.pp{jjj}(:,1);
                if (i > 1) || (jj > 1)
                    pp0{jj,i}{jjj}(:,1) = pp0{jj,i}{jjj}(:,1) +...
                        ((i - 1) * n_temporal + jj - 1) * params0_ud.Value(5);
                end
            end

            jj = jj + 1; % move to next saving spot

        end
        
        % record real-time ensemble data
        ensdata0.t(j+1,i) = min(pars0.delt) + ensdata0.t(j,i);
        ensdata0.n_agg(j+1,i) = length(pars0.n);
        ensdata0.tau(j+1,2*(i-1)+(1:2)) = [mean(pars0.tau), std(pars0.tau)];
        ensdata0.kn_kin(j+1,2*(i-1)+(1:2)) = [mean(pars0.kn_kin), std(pars0.kn_kin)];
        ensdata0.kn_diff(j+1,2*(i-1)+(1:2)) = [mean(pars0.kn_diff), std(pars0.kn_diff)];
        ensdata0.dpp(j+1,2*(i-1)+(1:2)) = [geomean(pars0.dpp_g(:,1)),...
            UTILS.GEOSTD(pars0.dpp_g(:,1))];
        ensdata0.sigmapp(j+1,2*(i-1)+(1:2)) = [geomean(pars0.dpp_g(:,2)),...
            UTILS.GEOSTD(pars0.dpp_g(:,2))];
        ensdata0.dm(j+1,2*(i-1)+(1:2)) = [geomean(pars0.dm),...
            UTILS.GEOSTD(pars0.dm)];
        if strcmp(opts0_mobil.mtd, 'interp')
            ensdata0.da(j+1,2*(i-1)+(1:2)) = [geomean(pars0.da),...
                UTILS.GEOSTD(pars0.da)];
        end

        UTILS.TEXTBAR([j, j_max]); % update progress textbar

        j = j + 1; % move to next iteration

    end

    % save primary particle data for the last step if there is...
    % ...storage space
    if jj <= n_temporal
        pp0{jj,i} = pars0.pp;
        pp0_n{jj,i} = cat(1, pars0.n);
        i_pp0{jj} = cell(length(cat(1, pars0.n)), 1);
        for jjj = 1 : length(cat(1, pars0.n))
            i_pp0{jj}{jjj} = pars0.pp{jjj}(:,1);
            if (i > 1) || (jj > 1)
                pp0{jj,i}{jjj}(:,1) = pp0{jj,i}{jjj}(:,1) +...
                    ((i - 1) * n_temporal + jj - 1) * params0_ud.Value(5);
            end
        end
    end

    % remove duplicate aggregates
    if jj > 1
        i_pp00 = cat(1, i_pp0{:});
        ij = nchoosek(1 : length(cell2mat(pp0_n(:,i))), 2);
        ind_rmv = [];
        for ii = 1 : length(ij)
            if isequal(mod(i_pp00{ij(ii,1)}(:,1), params0_ud.Value(5)),...
                    mod(i_pp00{ij(ii,2)}(:,1), params0_ud.Value(5)))
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

