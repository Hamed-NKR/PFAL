function [m_tls, b0_tls, ci_m, ci_b0] = SVD_LOGFIT(dist_log_2d)
% SVDFIT perfroms singular value decomposition to identify the spatial...
%   ...direction (if existing) of a given scattered set of data. Then...
%   ...a linear relation is chosen in log-log space corresponding...
%   ...to the smallest singular value. This relation yields the total...
%   ...least squares (TLS) fit. 95% confidence intervals are calculated...
%   ...based on the difference between the scattered data and TLS...
%   ...using a t-distribution.
% ----------------------------------------------------------------------- %
%
% Inputs:
%   dist_log_2d: Bivariate distribution in log-log space.
% ----------------------------------------------------------------------- %
%
% Outputs
%   m_tls: total least squares slope
%   b0_tls: total least squares intercept
%   ci_m: 95% confidence intervals for slope
%   ci_b0 95% confidence intervals for intercept
% ----------------------------------------------------------------------- %

[~, ~, V] = svd(dist_log_2d, 'econ'); % perform SVD on the data matrix
dir = V(:, end);  % direction vector with smallest singular value
m_tls = -dir(1) / dir(2);  % slope for total least square (TLS)
b0_tls = mean(dist_log_2d(:,2) - m_tls * dist_log_2d(:,1)); % TLS intercept

% calculate orthogonal residuals and variance
eps = (dist_log_2d(:,2) - (m_tls * dist_log_2d(:,1) + b0_tls)) ./...
    sqrt(1 + m_tls^2);
var_eps = var(eps);

% estimate standard errors for slope and intercept
n_tls = size(dist_log_2d, 1); % number of data points
mu_da = mean(dist_log_2d(:,1));
s_xx = sum((dist_log_2d(:,1) - mu_da).^2); % sum of squares of difference in x dir.
se_m = sqrt(var_eps / s_xx);  % standard error of slope
intercept_se = sqrt(var_eps * (1 / n_tls + mu_da^2 / s_xx)); % standard...
    % ...error of intercept

% calculate 95% confidence intervals
t_critical = tinv(0.975, n_tls - 2); % Critical t-value for 95% CI with...
    % ...(n-2) degrees of freedom
ci_m = [m_tls - t_critical * se_m, m_tls + t_critical * se_m]; % CI for slope
ci_b0 = [b0_tls - t_critical * intercept_se, b0_tls +...
    t_critical * intercept_se]; % CI for intercept

end

