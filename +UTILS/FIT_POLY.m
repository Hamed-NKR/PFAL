function [yfit, xfit, bounds_yfit, afit, bounds_afit, out] = FIT_POLY(x, y, group, resol, varargin)
% FIT_POLY fits a polynomial regression model (optionally weighted) and returns
% the fitted curve with uncertainty bounds, as well as the slope of the fitted
% curve in the predictor domain used for regression.
%
% The function is designed as a general interface so that additional fitting
% methods (e.g., ordinary least squares) can be added later via the 'Method'
% option. The current implementation supports a Bayesian linear model with a
% conjugate Normal–Inverse-Gamma prior.
% ----------------------------------------------------------------------- %
%
% Model form (in transformed space):
%   t = beta0 + beta1*z + beta2*z^2 + ... + betad*z^d + betag*group + eps
%
% where:
%   z = XTransform(x)   (default: log10)
%   t = YTransform(y)   (default: log10)
%
% and eps is Gaussian noise. When weights are provided, the fit is performed
% using importance weighting (weighted least squares likelihood), i.e. the
% objective emphasizes points with larger weights.
% ----------------------------------------------------------------------- %
%
% Inputs:
%   x       : independent variable in original (non-transformed) space
%   y       : dependent variable in original (non-transformed) space
%   group   : optional grouping covariate (same length as x); if empty/missing,
%             group = 0 (no grouping effect)
%   resol   : number of points used to evaluate the fitted curve; if empty/missing,
%             resol = 1000
%
% Name-value options (varargin):
%   'Degree'      : polynomial degree in transformed predictor z (default: 2)
%   'Weights'     : nonnegative importance weights (default: [] => all ones)
%   'Method'      : fitting backend (default: 'bayes'); placeholder for 'ols'
%   'XTransform'  : 'log10' or 'none' (default: 'log10')
%   'YTransform'  : 'log10' or 'none' (default: 'log10')
%   'NSamples'    : number of posterior draws for Bayesian backend (default: 1000)
%   'SlopeOffset' : constant added to the slope output (default: 0)
%
% Prior hyperparameters (Bayesian backend):
%   'PriorMu'     : prior mean of regression coefficients (default: zeros)
%   'PriorVScale' : prior covariance scale (default: 100, so V = 100*I)
%   'PriorA'      : inverse-gamma shape parameter for sigma^2 (default: 3)
%   'PriorB'      : inverse-gamma scale parameter for sigma^2 (default: 1)
% ----------------------------------------------------------------------- %
%
% Outputs:
%   xfit          : predictor samples in original x-units
%   yfit          : posterior mean of fitted response in original y-units
%   bounds_yfit   : 95% credible interval for fitted response in original y-units
%                  (two columns: [lower, upper])
%   afit          : posterior mean slope in transformed coordinates, evaluated at xfit,
%                  plus SlopeOffset
%   bounds_afit   : 95% credible interval for slope (two columns: [lower, upper]),
%                  plus SlopeOffset
%   out           : struct of diagnostic outputs (posterior params, draws, options)
% ----------------------------------------------------------------------- %

%% -------------------- input defaults and basic checks ------------------- %

% Handle missing optional positional inputs.
if nargin < 3 || isempty(group)
    group = zeros(size(x));
end
if nargin < 4 || isempty(resol)
    resol = 1000;
end

% Enforce column vectors.
x = x(:);
y = y(:);
group = group(:);

n = numel(x);
if numel(y) ~= n || numel(group) ~= n
    error('FIT_POLY:SizeMismatch', 'Inputs x, y, and group must have the same length.');
end

% ------------------------- parse name-value options ---------------------- %

opt = struct();
opt.Degree      = 2;
opt.Weights     = [];
opt.Method      = 'bayes';
opt.XTransform  = 'log10';
opt.YTransform  = 'log10';
opt.NSamples    = 1000;
opt.SlopeOffset = 0;

% Prior defaults match the original function spirit: diffuse Normal prior and
% moderately informative inverse-gamma on sigma^2.
opt.PriorMu     = [];    % if empty, set to zeros with correct length later
opt.PriorVScale = 100;
opt.PriorA      = 3;
opt.PriorB      = 1;

if ~isempty(varargin)
    if mod(numel(varargin),2) ~= 0
        error('FIT_POLY:BadVarargin', 'Name-value inputs must come in pairs.');
    end
    for k = 1:2:numel(varargin)
        name = varargin{k};
        value = varargin{k+1};
        if ~ischar(name) && ~isstring(name)
            error('FIT_POLY:BadOptionName', 'Option names must be strings.');
        end
        name = char(lower(string(name)));

        switch name
            case 'degree'
                opt.Degree = value;
            case 'weights'
                opt.Weights = value;
            case 'method'
                opt.Method = char(lower(string(value)));
            case 'xtransform'
                opt.XTransform = char(lower(string(value)));
            case 'ytransform'
                opt.YTransform = char(lower(string(value)));
            case 'nsamples'
                opt.NSamples = value;
            case 'slopeoffset'
                opt.SlopeOffset = value;

            case 'priormu'
                opt.PriorMu = value;
            case 'priorvscale'
                opt.PriorVScale = value;
            case 'priora'
                opt.PriorA = value;
            case 'priorb'
                opt.PriorB = value;

            otherwise
                error('FIT_POLY:UnknownOption', 'Unknown option: %s', name);
        end
    end
end

% Validate key options.
deg = opt.Degree;
if ~isscalar(deg) || deg < 1 || deg ~= floor(deg)
    error('FIT_POLY:BadDegree', 'Degree must be a positive integer (>=1).');
end

if ~isscalar(resol) || resol < 2 || resol ~= floor(resol)
    error('FIT_POLY:BadResol', 'resol must be an integer >= 2.');
end

nsamp = opt.NSamples;
if ~isscalar(nsamp) || nsamp < 100 || nsamp ~= floor(nsamp)
    error('FIT_POLY:BadNSamples', 'NSamples must be an integer >= 100.');
end

%% ----------------------------- transforms -------------------------------- %

% Define forward and inverse transforms for x and y.
% The regression always occurs in the transformed space (z, t).
switch opt.XTransform
    case 'log10'
        if any(x <= 0)
            error('FIT_POLY:NonPositiveX', 'x must be > 0 when XTransform = log10.');
        end
        fX  = @(u) log10(u);
        iX  = @(u) 10.^u;
    case 'none'
        fX  = @(u) u;
        iX  = @(u) u;
    otherwise
        error('FIT_POLY:BadXTransform', 'XTransform must be ''log10'' or ''none''.');
end

switch opt.YTransform
    case 'log10'
        if any(y <= 0)
            error('FIT_POLY:NonPositiveY', 'y must be > 0 when YTransform = log10.');
        end
        fY  = @(u) log10(u);
        iY  = @(u) 10.^u;
    case 'none'
        fY  = @(u) u;
        iY  = @(u) u;
    otherwise
        error('FIT_POLY:BadYTransform', 'YTransform must be ''log10'' or ''none''.');
end

% Apply transforms.
z = fX(x);     % transformed predictor
t = fY(y);     % transformed response

% Warn when a nonzero slope offset is used outside log–log context.
if opt.SlopeOffset ~= 0
    if ~(strcmp(opt.XTransform,'log10') && strcmp(opt.YTransform,'log10'))
        warning('FIT_POLY:SlopeOffsetNonLog', ...
            ['SlopeOffset is applied to slopes computed in transformed coordinates. ' ...
             'Interpretation differs from a log–log slope unless both transforms are log10.']);
    end
end

%% ----------------------------- weights ----------------------------------- %

% Importance weights enter through the likelihood as a weighted sum of squared
% residuals. This is equivalent to multiplying each row of the design matrix
% and response by sqrt(w_i).
if isempty(opt.Weights)
    w = ones(n,1);
else
    w = opt.Weights(:);
    if numel(w) ~= n
        error('FIT_POLY:BadWeights', 'Weights must have the same length as x and y.');
    end
    if any(~isfinite(w)) || any(w <= 0)
        error('FIT_POLY:BadWeights', 'Weights must be finite and strictly positive.');
    end
end

% Normalize weights by their mean for numerical conditioning. This preserves
% relative importance while avoiding unnecessary scaling of the noise variance.
w = w ./ mean(w);

sqrtw = sqrt(w);

%% ---------------------- build design matrices ---------------------------- %

% Design matrix in transformed space:
%   Xfull = [1, z, z^2, ..., z^deg, group]
%
% The final column includes the group covariate, consistent with the earlier
% implementation that used the group mean for prediction.
Xpoly = zeros(n, deg);
for k = 1:deg
    Xpoly(:,k) = z.^k;
end

Xfull = [ones(n,1), Xpoly, group];   % size: n x K
K = size(Xfull,2);

% Prior mean and covariance for beta (includes intercept).
if isempty(opt.PriorMu)
    Mu = zeros(K,1);
else
    Mu = opt.PriorMu(:);
    if numel(Mu) ~= K
        error('FIT_POLY:BadPriorMu', 'PriorMu must have length %d (intercept + degree terms + group).', K);
    end
end

Vscale = opt.PriorVScale;
if ~isscalar(Vscale) || Vscale <= 0
    error('FIT_POLY:BadPriorVScale', 'PriorVScale must be a positive scalar.');
end
V = Vscale * eye(K);

A0 = opt.PriorA;
B0 = opt.PriorB;
if ~isscalar(A0) || A0 <= 0 || ~isscalar(B0) || B0 <= 0
    error('FIT_POLY:BadPriorAB', 'PriorA and PriorB must be positive scalars.');
end

%% ----------------------------- method ------------------------------------ %

switch opt.Method
    case 'bayes'
        % Bayesian conjugate regression with Normal–Inverse-Gamma prior.
        % The weighted likelihood is implemented via row scaling by sqrt(w).

    case 'ols'
        % Placeholder for future expansion. The design/prediction plumbing in
        % this function is organized so that OLS can be added with minimal
        % duplication (beta estimate + covariance -> confidence intervals).
        error('FIT_POLY:MethodNotImplemented', ...
              'Method ''ols'' is not implemented yet. Use Method ''bayes'' for now.');

    otherwise
        error('FIT_POLY:BadMethod', 'Unknown Method: %s', opt.Method);
end

%% ------------------ Bayesian posterior (closed form) ---------------------- %

% Weighted regression via row scaling:
Xw = Xfull .* sqrtw;   % each row i multiplied by sqrt(w_i)
tw = t .* sqrtw;

% Posterior parameters for Normal–Inverse-Gamma model:
% Prior:
%   beta | sigma^2 ~ N(Mu, sigma^2 * V)
%   sigma^2 ~ InvGamma(A0, B0)
%
% Likelihood (weighted):
%   tw ~ N(Xw*beta, sigma^2 I)
%
% Posterior:
%   beta | sigma^2, data ~ N(Mun, sigma^2 * Vn)
%   sigma^2 | data ~ InvGamma(An, Bn)

invV = inv(V);

Vn = inv(invV + (Xw' * Xw));
Mun = Vn * (invV * Mu + (Xw' * tw));

An = A0 + n/2;

% Standard conjugate update for Bn:
% Bn = B0 + 0.5 * ( t'W t + Mu'V^{-1}Mu - Mun'Vn^{-1}Mun )
%
% where t'W t is implemented as tw'*tw.
invVn = inv(Vn);
Bn = B0 + 0.5 * ( (tw' * tw) + (Mu' * invV * Mu) - (Mun' * invVn * Mun) );

%% -------------------------- posterior sampling --------------------------- %

% Sample sigma^2 from InvGamma(An, Bn).
% If sigma^2 ~ InvGamma(a, b) in the "scale" parameterization:
%   p(sigma^2) ∝ (sigma^2)^(-a-1) exp(-b/sigma^2)
% then tau = 1/sigma^2 ~ Gamma(a, scale = 1/b).
tau = gamrnd(An, 1./Bn, [1, nsamp]);
sigma2_draws = 1 ./ tau;

% Sample beta | sigma^2 ~ N(Mun, sigma^2 * Vn).
% Use Cholesky factorization for stable sampling.
L = chol(Vn, 'lower');                 % Vn = L*L'
Zrand = randn(K, nsamp);               % standard normal draws
beta_draws = Mun + (L * Zrand) .* sqrt(sigma2_draws);  % broadcast sqrt(sigma2) across rows

%% ----------------------- prediction grid & fit --------------------------- %

% Construct evaluation grid in transformed predictor space.
zfit = linspace(min(z), max(z), resol)';

% Set group to its mean for fitted curve evaluation (matches earlier behavior).
gbar = mean(group);

Xpoly_fit = zeros(resol, deg);
for k = 1:deg
    Xpoly_fit(:,k) = zfit.^k;
end

Xfit_full = [ones(resol,1), Xpoly_fit, gbar * ones(resol,1)];   % resol x K

% Posterior draws of fitted mean response in transformed space:
% t_pred_samples(:,s) = Xfit_full * beta_draws(:,s)
t_pred_samples = Xfit_full * beta_draws;

% Summarize fitted curve in transformed space.
tfit_mean  = mean(t_pred_samples, 2);
tfit_lower = prctile(t_pred_samples, 2.5, 2);
tfit_upper = prctile(t_pred_samples, 97.5, 2);

% Transform back to original space.
xfit = iX(zfit);
yfit = iY(tfit_mean);
bounds_yfit = [iY(tfit_lower), iY(tfit_upper)];

%% ---------------------- slope (derivative) output ------------------------ %

% The slope is computed as dt/dz for the fitted polynomial in transformed space:
%   dt/dz = sum_{k=1..deg} k * beta_k * z^(k-1)
%
% Coefficients beta_1..beta_deg correspond to rows 2..(deg+1) of beta_draws.
beta_poly_draws = beta_draws(2:(deg+1), :);    % deg x nsamp

Dbasis = zeros(resol, deg);
for k = 1:deg
    Dbasis(:,k) = k * (zfit.^(k-1));
end

% Derivative draws in transformed coordinates.
slope_samples = Dbasis * beta_poly_draws;      % resol x nsamp

% Summarize slope and apply slope offset.
afit = mean(slope_samples, 2) + opt.SlopeOffset;
bounds_afit = [prctile(slope_samples, 2.5, 2), prctile(slope_samples, 97.5, 2)] + opt.SlopeOffset;

% Informational warning when not log–log. The slope is still returned, but its
% meaning is the derivative in the transformed coordinates (t vs z), not
% necessarily a log–log exponent.
if ~(strcmp(opt.XTransform,'log10') && strcmp(opt.YTransform,'log10'))
    warning('FIT_POLY:SlopeInterpretation', ...
        ['Slope is reported as dt/dz in transformed coordinates. ' ...
         'Interpretation depends on XTransform and YTransform.']);
end

%% ---------------------------- diagnostics -------------------------------- %

out = struct();
out.Options = opt;

out.Transforms = struct();
out.Transforms.z = z;
out.Transforms.t = t;
out.Transforms.zfit = zfit;

out.Design = struct();
out.Design.Xfull = Xfull;
out.Design.Xfit_full = Xfit_full;

out.Posterior = struct();
out.Posterior.Mu  = Mu;
out.Posterior.V   = V;
out.Posterior.A0  = A0;
out.Posterior.B0  = B0;
out.Posterior.Mun = Mun;
out.Posterior.Vn  = Vn;
out.Posterior.An  = An;
out.Posterior.Bn  = Bn;

out.Draws = struct();
out.Draws.beta = beta_draws;
out.Draws.sigma2 = sigma2_draws;

out.Weights = w;

end
