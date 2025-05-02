function [dm_fit, rho_fit_mean, rho_fit_lower, rho_fit_upper] = OLSLOGPOLYFIT2(dm_data, rho_data)
% OLSLOGPOLYFIT2: Log-log 2nd order OLS regression with 95% CI (CI is buggy)
% ----------------------------------------------------------------------- %
% Inputs: dm_data [m], rho_data [kg/m^3]
% ----------------------------------------------------------------------- %
% Outputs: fitted diameter and CI curves [in meters and kg/m^3]
% ----------------------------------------------------------------------- %

    % Filter and log-transform
    valid = isfinite(dm_data) & isfinite(rho_data) & dm_data > 0 & rho_data > 0;
    x_log = log10(dm_data(valid));
    y_log = log10(rho_data(valid));

    % Center x for numerical stability
    x0 = mean(x_log);
    x_c = x_log - x0;

    % Fit OLS with centered predictor
    X = [ones(size(x_c)), x_c, x_c.^2];
    beta = X \ y_log;

    % Residual stats
    n = length(y_log);
    p = size(X,2);
    res = y_log - X * beta;
    sigma2 = sum(res.^2) / (n - p);
    C = inv(X' * X);
    tval = tinv(0.975, n - p);

    % Generate predictions on same centered domain
    x_fit_log = linspace(min(x_log), max(x_log), 500)';
    x_fit_c = x_fit_log - x0;
    X_fit = [ones(size(x_fit_c)), x_fit_c, x_fit_c.^2];

    % Prediction and CI in log-space
    y_pred = X_fit * beta;
    se = sqrt(sum((X_fit * C) .* X_fit, 2) * sigma2);
    y_lower = y_pred - tval * se;
    y_upper = y_pred + tval * se;

    % Convert to linear domain
    dm_fit = 10.^x_fit_log;
    rho_fit_mean = 10.^y_pred;
    rho_fit_lower = 10.^y_lower;
    rho_fit_upper = 10.^y_upper;

end
