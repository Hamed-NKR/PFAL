function [yfit, xfit, bounds_yfit, afit, bounds_afit] =...
    BAYESFIT_POLY4(x, y, group, resol)
% BAYESFIT perfroms a fourth-order polynomial regression in log-log...
%   ...space based on bayesian inference. This is primarily used for...
%   ...capturing curvature in the trends of effective density vs. ...
%   ...mobility diameter. Means and 95% confidence intervals are output...
%   ...from this function for both the main data in original space and...
%   ...slope in log space. Effect of grouping on the data is considered.
% ----------------------------------------------------------------------- %
%
% Inputs:
%   x: independant variable in original (non-logarithmic) space
%   y: depend variable (non-log) with dispersion (to be regressed on x)
%   group: an indentifier for different groups existing within x
%   resol: resolution of output fit
% ----------------------------------------------------------------------- %
%
% Outputs
%   xfit: predictor samples
%   yfit: mean of response samples
%   bounds_yfit: 95% confidence intervals for the response samples(a...
%       ...two-column array with the first column showing lower bounds...
%       ...and second one showing upper bounds)
%   afit: mean of slopes of respose sample
%   bounds_afit: 95% confidence intervals for the slopes of response...
%       ...samples(again: [lower bounds, upper bounds])
% ----------------------------------------------------------------------- %

%% Bayesian polynomial regression

% prepare predictors and response
log_x = log10(x);
log_y = log10(y);

X = [log_x, log_x.^2, log_x.^3, log_x.^4, group];  % a fourth-order...
    % ...polynomial ensures better capturing of curvature

% define prior
p = size(X,2); % number of predictors (excluding intercept)
Mu = zeros(p + 1, 1); % prior mean (intercept + p)
V = 100 * eye(p + 1);  % prior covariance
A = 3; B = 1; % inverse gamma prior on sigma^2

% create semiconjugate prior model
Mdl = semiconjugateblm(p, 'Intercept', true, 'Mu', Mu, 'V', V, 'A', A,...
    'B', B);

% estimate posterior
PosteriorMdl = estimate(Mdl, X, log_y);  % uses Gibbs sampling

% get posterior samples
B_samples = PosteriorMdl.BetaDraws;

% prediction setup
log_xfit = linspace(min(log_x), max(log_x), resol)';
mu_group = mean(group);
X_fit = [ones(size(log_xfit)), log_xfit, log_xfit.^2, log_xfit.^3, log_xfit.^4,...
    mu_group * ones(size(log_xfit))];

% predict y from all posterior samples
y_pred_samples = X_fit * B_samples;

% compute posterior predictive summary
log_yfit = mean(y_pred_samples, 2);
y_lower = prctile(y_pred_samples, 2.5, 2);
y_upper = prctile(y_pred_samples, 97.5, 2);

% transform back to linear space
xfit = 10.^log_xfit;
yfit = 10.^log_yfit;
bounds_yfit = [10.^y_lower, 10.^y_upper];

%% calculate slope from regression results (derivative of polynomial)

% allocate derivative (slope)
n_points = length(log_xfit);
n_samples = size(B_samples, 2);
alpha_samples = zeros(n_points, n_samples);

% extract polynomial coefficients (skipping intercept and group term)
b1 = B_samples(2, :);
b2 = B_samples(3, :);
b3 = B_samples(4, :);
b4 = B_samples(5, :);

% compute polynomial slopes at each x_fit using all posterior samples
for i = 1 : n_samples
    alpha_samples(:, i) = b1(i) + 2 * b2(i) * log_xfit +...
        3 * b3(i) * log_xfit.^2 + 4 * b4(i) * log_xfit.^3;
end

% posterior summary of slope
afit = mean(alpha_samples, 2);
bounds_afit = [prctile(alpha_samples, 2.5, 2),...
    prctile(alpha_samples, 97.5, 2)];

end

