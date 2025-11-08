clear; clc; close all;

load('D:\Hamed\CND\PhD\Weekly\2025\02OCT2025\Realizations_02OCT25\tem_manual.mat');

% fdir_valid = 'D:\Hamed\CND\PhD\Weekly\2025\02OCT2025\Realizations_02OCT25\real4\Valid_02-Oct-2025_13-08-51';
% fname_valid = 'Valid_02-Oct-2025_13-08-51.mat';
fdir_valid = 'D:\Hamed\CND\PhD\Weekly\2025\02OCT2025\Realizations_02OCT25\real3\Valid_02-Oct-2025_12-05-58';
fname_valid = 'Valid_02-Oct-2025_12-05-58.mat';
vars_valid = {'bayesfit_valid', 'dm_uc', 'rho_eff_uc', 'f5'};
load(fullfile(fdir_valid, fname_valid), vars_valid{:});

hold on

rho_eff_convert = TRANSP.CONVERT_DPP([ci_ens_lal_unw(1),...
    gm_ens_lal_unw, ci_ens_lal_unw(2)], 'table');

manualShade_x = [ci_da_lal_unw(1) ci_da_lal_unw(2)...
    ci_da_lal_unw(2) ci_da_lal_unw(1)];
manualShade_y = [min(rho_eff_convert) min(rho_eff_convert)...
    max(rho_eff_convert) max(rho_eff_convert)];

fill(manualShade_x, manualShade_y, 'r', 'FaceAlpha', 0.3,...
    'FaceColor', 'k', 'EdgeColor', 'k');

lgds = findobj(f5, 'Type', 'Legend');
lgds.String{6} = 'Manual TEM data converted';

exportgraphics(f5, strcat('outputs\', 'simul-vs-exp-v3.jpg'),...
    'BackgroundColor','none', 'Resolution', 300)

