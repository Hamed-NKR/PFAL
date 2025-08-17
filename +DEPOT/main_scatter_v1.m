clc
clear
close all
warning('off')

f1 = figure;
f1.Position = [100, 100, 600, 600];
set(f1, 'color', 'white');

% plot universal correlation
r0 = (2e4 / 1e0) ^ (1 / (1e4 - 1));
da0 = 1e0 * ones(1e4, 1) .* r0 .^ (((1:1e4)-1)');
dpp0 = 17.8 * (da0 / 100) .^ (0.35);
plt_0 = plot(da0, dpp0, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 2);
hold on

load('D:\Hamed\CND\PhD\My Articles\DLCA1\Results\DAT\database.mat', 'parsdata_sigma')

da_uc = 1e9 * parsdata_sigma{4}(1).da;
dpp_uc = 1e9 * parsdata_sigma{4}(1).dpp_g(:,1);
plt_uc = scatter(da_uc, dpp_uc, 15, hex2rgb('#295F98'), 'o', 'LineWidth', 1);

% Generate distances that follow a normal distribution (mean 0, std deviation 1)
% mu_g = 0;
% sigma_g = 2;
% mu = mu_g * exp((sigma_g^2)/2);
% sigma = mu^2 * (exp(sigma_g^2) - 1);
mu = 0.3;
sigma = 0.1;
dist = normrnd(mu, sigma, size(da_uc));

% Compute the perpendicular direction vector to the line
m = 0.35;
perp_da = -1 / sqrt(1 + m^2);  % x-component of the perpendicular direction
perp_dpp = m / sqrt(1 + m^2);   % y-component of the perpendicular direction

% Scatter the points perpendicularly from the line using the normal distances
da_scat = exp(log(da_uc) + dist * perp_da);
dpp_scat = exp(log(dpp_uc) + dist * perp_dpp);

plt_scat = scatter(da_scat, dpp_scat, 15, hex2rgb('#EF9C66'), '^', 'LineWidth', 1);

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([10,1500])
ylim([5,50])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
legend(cat(2, plt_0, plt_uc, plt_scat),...
    cat(2, {'Olfert and Rogak (2019)'},...
    {'Original scaling'}, {'Secondary scattering'}),...
    'interpreter', 'latex', 'FontSize', 12, 'location', 'northwest')

