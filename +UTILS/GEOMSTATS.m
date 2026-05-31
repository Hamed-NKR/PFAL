function [gmean, gstd, ci95] = GEOMSTATS(x, w)
%GEOMSTATS Compute geometric mean, geometric std. dev., and 95% CI.
%   [GMEAN, GSTD, CI95] = UTILS.GEOMSTATS(X) computes unweighted
%   log-normal statistics for positive finite values in X.
%
%   [GMEAN, GSTD, CI95] = UTILS.GEOMSTATS(X, W) computes weighted
%   log-normal statistics using positive finite weights W and Kish
%   effective sample size for the confidence interval.
%
%   Adapted from Hamed-NKR/atems +morph/geomstats.m on branch HN.

x = x(:);

if nargin < 2 || isempty(w)
    x = x(isfinite(x) & x > 0);
    n = numel(x);

    if n == 0
        gmean = NaN;
        gstd = NaN;
        ci95 = [NaN NaN];
        return
    elseif n == 1
        gmean = x;
        gstd = 1;
        ci95 = [x x];
        return
    end

    y = log(x);
    mu = mean(y);
    s = std(y, 0);

    gmean = exp(mu);
    gstd = exp(s);

    se = s / sqrt(n);
    z = 1.96;
    ci95 = exp(mu + z * [-1 1] * se);
else
    w = w(:);
    ok = isfinite(x) & x > 0 & isfinite(w) & w > 0;
    x = x(ok);
    w = w(ok);

    if isempty(x)
        gmean = NaN;
        gstd = NaN;
        ci95 = [NaN NaN];
        return
    end

    w = w / sum(w);
    y = log(x);
    mu = sum(w .* y);

    gmean = exp(mu);
    variance_log = sum(w .* (y - mu) .^ 2);
    gstd = exp(sqrt(variance_log));

    n_eff = 1 / sum(w .^ 2);
    if n_eff > 1
        se = sqrt(variance_log) / sqrt(n_eff);
        z = 1.96;
        ci95 = exp(mu + z * [-1 1] * se);
    else
        ci95 = [gmean gmean];
        warning('PFAL:GEOMSTATS:LowEffectiveN', ...
            'Effective sample size n_eff = %.2f <= 1. CI is not reliable.', ...
            n_eff);
    end
end

end
