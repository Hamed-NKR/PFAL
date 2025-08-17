%%% A script to render individual aggregates along with their...
%%% ...characteristic sizes

clc
clear
clf('reset')
close all
warning('off')

%% initialize %%

% address of data library to be imported
fdir = 'D:\Hamed\CND\PhD\My Articles\DLCA2\mainscatter_sigmapp13\SCAT';
fname = 'SCAT-18NOV24';
varname = 'pars_out';

ind_agg = [];

rho0 = 1860; % material density

%% load and prepare the data %%

% load a pars stcructure containg aggregate data
load(strcat(fdir, '\', fname, '.mat'), varname)

eval(['pars_render' ' = ' varname ';']); % rename loaded structure

eval(['clear ', varname]) % remove older name

if ~isfield(pars_render, 'pp')
    disp(' ')
    error('Library does not contain aggregates!')
end

% pars_render = rmfield(pars_render, {'n', 'dpp_g', 'da', 'dg', 'dv'});

n_agg = size(pars_render.pp,1); % number of aggregates in the pars structure

% calculate number of primary particles within each aggregate if not...
% ...existing
if ~isfield(pars_render, 'n')
    pars_render.n = zeros(n_agg,1);

    for i = 1 : n_agg
        pars_render.n(i) = size(pars_render.pp{i}, 1);
    end
end

% calculate characteristic sizes if not available
if ~isfield(pars_render, 'dpp_g') || ~isfield(pars_render, 'da') ||...
        ~isfield(pars_render, 'dv') || ~isfield(pars_render, 'dg')
    pars_render = PAR.SIZING(pars_render);
    pars_render.da = 2 * sqrt(PAR.PROJECTION(pars_render, [], 1e3, 10) / pi);
end

% calculate mobility diameter if needed
[params_ud, params_const] = TRANSP.INIT_PARAMS('LD2_Params');
[~, fl] = TRANSP.INIT_DOM(params_ud, params_const);
pars_render.dm = TRANSP.DIAMOBIL(pars_render.dg, pars_render.da, fl);

% calculate aerodynamic diameter
pars_render.dae = TRANSP.AERODYN(pars_render.pp, pars_render.dm, rho0);

%% visualize aggregate strcuture %%

f0 = figure(7);
f0.Position = [0, 0, 600, 600];
set(f0, 'color', 'white');

% select a random aggregate if not already selected
if isempty(ind_agg)
    ind_agg = randperm(n_agg, 1);
end

% adjust aggregate color
opts_rnd.cm = sky;
% opts_rnd.cm = flip(opts_rnd.cm);
opts_rnd.cloc = 0.4;

% render the aggregate
plt0_agg = UTILS.PLOTPP(pars_render.pp{ind_agg}(:,3),...
    pars_render.pp{ind_agg}(:,4), pars_render.pp{ind_agg}(:,5),...
    pars_render.pp{ind_agg}(:,2), [], opts_rnd);


%% draw primary particle diameter %%

% determine center of mass
r_com = PAR.COM(pars_render.pp(ind_agg), pars_render.n(ind_agg));

f1 = figure(1);
f1.Position = [50, 50, 600, 600];
set(f1, 'color', 'white');

% select a random aggregate if not already selected
if isempty(ind_agg)
    ind_agg = randperm(n_agg, 1);
end

% render the aggregate
plt1_agg = UTILS.PLOTPP(pars_render.pp{ind_agg}(:,3),...
    pars_render.pp{ind_agg}(:,4), pars_render.pp{ind_agg}(:,5),...
    pars_render.pp{ind_agg}(:,2), [], opts_rnd);

hold on

% get geometric mean and sd
dpp = geomean(pars_render.dpp_g);
sigmapp = dpp(2);
dpp = dpp(1);

[X,Y,Z] = sphere(60); % baseline sphere to draw characteristic sizes

% draw primary particle size
plt_dpp = surf(X * dpp / 2 + r_com(1), Y .* dpp / 2 + r_com(2),...
    Z .* dpp / 2 + r_com(3), 'EdgeColor', 'none', 'FaceColor',...
    hex2rgb('#8174A0'), 'FaceAlpha', 0.3);
plt_dpp = setgraphics(plt_dpp);

% draw upper extent based on GSD
plt_dpp_max = surf(X * (dpp * sigmapp) / 2 + r_com(1), Y .* (dpp * sigmapp)...
    / 2 + r_com(2), Z .* (dpp * sigmapp) / 2 + r_com(3), 'EdgeColor', 'none',...
    'FaceColor', hex2rgb('#8174A0'), 'FaceAlpha', 0.3);
plt_dpp_max = setgraphics(plt_dpp_max);

% draw lower extent based on GSD
plt_dpp_min = surf(X * (dpp / sigmapp) / 2 + r_com(1), Y .* (dpp / sigmapp)...
    / 2 + r_com(2), Z .* (dpp / sigmapp) / 2 + r_com(3), 'EdgeColor', 'none',...
    'FaceColor', hex2rgb('#8174A0'), 'FaceAlpha', 0.3);
plt_dpp_min = setgraphics(plt_dpp_min);

%% draw projected area diameter %%

f2 = figure(2);
f2.Position = [100, 100, 600, 600];
set(f2, 'color', 'white');

plt2_agg = UTILS.PLOTPP(pars_render.pp{ind_agg}(:,3),...
    pars_render.pp{ind_agg}(:,4), pars_render.pp{ind_agg}(:,5),...
    pars_render.pp{ind_agg}(:,2), [], opts_rnd);

hold on

plt_da = surf(X * pars_render.da(ind_agg) / 2 + r_com(1), Y .*...
    pars_render.da(ind_agg) / 2 + r_com(2), Z .* pars_render.da(ind_agg) /...
    2 + r_com(3), 'EdgeColor', 'none', 'FaceColor', hex2rgb('#8174A0'),...
    'FaceAlpha', 0.3);
plt_da = setgraphics(plt_da);

%% draw gyration diameter %%

f3 = figure(3);
f3.Position = [150, 150, 600, 600];
set(f3, 'color', 'white');

plt3_agg = UTILS.PLOTPP(pars_render.pp{ind_agg}(:,3),...
    pars_render.pp{ind_agg}(:,4), pars_render.pp{ind_agg}(:,5),...
    pars_render.pp{ind_agg}(:,2), [], opts_rnd);

hold on

plt_dg = surf(X * pars_render.dg(ind_agg) / 2 + r_com(1), Y .*...
    pars_render.dg(ind_agg) / 2 + r_com(2), Z .* pars_render.dg(ind_agg) /...
    2 + r_com(3), 'EdgeColor', 'none', 'FaceColor', hex2rgb('#8174A0'),...
    'FaceAlpha', 0.3);
plt_dg = setgraphics(plt_dg);

%% draw volume diameter %%

f4 = figure(4);
f4.Position = [200, 200, 600, 600];
set(f4, 'color', 'white');

plt4_agg = UTILS.PLOTPP(pars_render.pp{ind_agg}(:,3),...
    pars_render.pp{ind_agg}(:,4), pars_render.pp{ind_agg}(:,5),...
    pars_render.pp{ind_agg}(:,2), [], opts_rnd);

hold on

plt_dv = surf(X * pars_render.dv(ind_agg) / 2 + r_com(1), Y .*...
    pars_render.dv(ind_agg) / 2 + r_com(2), Z .* pars_render.dv(ind_agg) /...
    2 + r_com(3), 'EdgeColor', 'none', 'FaceColor', hex2rgb('#8174A0'),...
    'FaceAlpha', 0.3);
plt_dv = setgraphics(plt_dv);

%% draw mobility diameter %%

f5 = figure(5);
f5.Position = [250, 250, 600, 600];
set(f5, 'color', 'white');

plt5_agg = UTILS.PLOTPP(pars_render.pp{ind_agg}(:,3),...
    pars_render.pp{ind_agg}(:,4), pars_render.pp{ind_agg}(:,5),...
    pars_render.pp{ind_agg}(:,2), [], opts_rnd);

hold on

plt_dm = surf(X * pars_render.dm(ind_agg) / 2 + r_com(1), Y .*...
    pars_render.dm(ind_agg) / 2 + r_com(2), Z .* pars_render.dm(ind_agg) /...
    2 + r_com(3), 'EdgeColor', 'none', 'FaceColor', hex2rgb('#8174A0'),...
    'FaceAlpha', 0.3);
plt_dm = setgraphics(plt_dm);

%% draw aerodynamic diameter %%

f6 = figure(6);
f6.Position = [300, 300, 600, 600];
set(f6, 'color', 'white');

plt6_agg = UTILS.PLOTPP(pars_render.pp{ind_agg}(:,3),...
    pars_render.pp{ind_agg}(:,4), pars_render.pp{ind_agg}(:,5),...
    pars_render.pp{ind_agg}(:,2), [], opts_rnd);

hold on

plt_dae = surf(X * pars_render.dae(ind_agg) / 2 + r_com(1), Y .*...
    pars_render.dae(ind_agg) / 2 + r_com(2), Z .* pars_render.dae(ind_agg) /...
    2 + r_com(3), 'EdgeColor', 'none', 'FaceColor', hex2rgb('#8174A0'),...
    'FaceAlpha', 0.3);
plt_dae = setgraphics(plt_dae);

%% Export the rendering %%

if ~isfolder('outputs')
    mkdir('outputs'); % If the directory doesn't exist, create it 
end

exportgraphics(f0, 'outputs\render0.png', 'BackgroundColor',...
    'none', 'ContentType', 'vector', 'Resolution', 300) % save dpp figure

exportgraphics(f1, 'outputs\render_dpp.png', 'BackgroundColor',...
    'none', 'ContentType', 'vector', 'Resolution', 300) % save dpp figure

exportgraphics(f2, 'outputs\render_da.png', 'BackgroundColor',...
    'none', 'ContentType', 'vector', 'Resolution', 300) % save da figure

exportgraphics(f3, 'outputs\render_dg.png', 'BackgroundColor',...
    'none', 'ContentType', 'vector', 'Resolution', 300) % save dg figure

exportgraphics(f4, 'outputs\render_dv.png', 'BackgroundColor',...
    'none', 'ContentType', 'vector', 'Resolution', 300) % save dv figure

exportgraphics(f5, 'outputs\render_dm.png', 'BackgroundColor',...
    'none', 'ContentType', 'vector', 'Resolution', 300) % save dm figure

exportgraphics(f6, 'outputs\render_dae.png', 'BackgroundColor',...
    'none', 'ContentType', 'vector', 'Resolution', 300) % save dae figure


%% on-deman functionality %%

function plt = setgraphics(plt)

    plt.FaceLighting = 'gouraud';
    plt.AmbientStrength = 0.8;
    plt.DiffuseStrength = 0.2;
    plt.SpecularStrength = 0.05;
    plt.SpecularExponent = 2;
    plt.BackFaceLighting = 'lit';

end