function [rho_eff, dm_nm] = CONVERT_DPP(dp_nm, method, params)
% "CONVERT_DPP" maps primary particle diameter to effective density and...
%   ...mobility diameter based on correlations/data discussed in...
%   ...Sipkens et al. (2023): "Overview of methods...".
% ----------------------------------------------------------------------- %
%   [rho_eff, dm_nm] = dp_to_rhoeff_dm(dp_nm)
%   [rho_eff, dm_nm] = dp_to_rhoeff_dm(dp_nm, method)
%   [rho_eff, dm_nm] = dp_to_rhoeff_dm(dp_nm, method, params)
% ----------------------------------------------------------------------- %
% Inputs
%   dp_nm  : Nx1 primary particle diameters [nm]
%   method : 'equations' (default) or 'table'
%   params : struct to override defaults (fields depend on method)
% ----------------------------------------------------------------------- %
% Outputs
%   rho_eff : Nx1 effective density [kg/m^3]
%   dm_nm   : Nx1 mobility diameter [nm]
% -------------------------------------------------------------------------
% PHYSICS USED (equations method)
% 1) TEM scaling:        dp = dp100 * (dm/100)^D_TEM   (using dm ≈ d_A)
%    => dm = 100 * (dp/dp100)^(1/D_TEM)
% 2) Mass–mobility law:  m = m100 * (dm/100)^D_m
% 3) Effective density:  rho_eff = 6*m / (pi*dm^3)     (SI units)
% -------------------------------------------------------------------------
% Defaults (typical fresh soot in air):
%   D_m    = 2.48
%   m100   = 0.267 fg (mass at dm = 100 nm)
%   D_TEM  = 0.35
%   dp100  = 17.8 nm  (primary size at dm = 100 nm)
% -------------------------------------------------------------------------
% Usage:
%   Example primary sizes:
%   dp = [12; 18; 25; 35];   % nm
% 
%   1) Using equations (smooth, parametric; recommended)
%   [rho1, dm1] = dp_to_rhoeff_dm(dp, 'equations');
% 
%   2) Using the table (discrete, matches table 2 in Sipkens et al. (2023))
%   [rho2, dm2] = dp_to_rhoeff_dm(dp, 'table');
% -------------------------------------------------------------------------

if nargin < 2 || isempty(method), method = 'equations'; end
if nargin < 3, params = struct(); end

dp_nm = dp_nm(:);                   % ensure column

switch lower(method)
    case 'equations'
        % Defaults (override with params if provided)
        D_m   = getf(params, 'D_m',   2.48);
        m100_fg = getf(params, 'm100_fg', 0.267);  % fg
        D_TEM = getf(params, 'D_TEM', 0.35);
        dp100 = getf(params, 'dp100', 17.8);       % nm

        % 1) dp -> dm using inverted TEM scaling (dm ≈ d_A)
        dm_nm = 100 .* (dp_nm ./ dp100) .^ (1./D_TEM);

        % 2) mass from mass–mobility (keep nm & fg here for stability)
        m_fg  = m100_fg .* (dm_nm ./ 100) .^ D_m;  % fg

        % 3) effective density in SI
        dm_m  = dm_nm * 1e-9;                      % m
        m_kg  = m_fg  * 1e-15;                     % kg
        rho_eff = (6 .* m_kg) ./ (pi .* dm_m.^3);  % kg/m^3

    case 'table'
        % Digitized values from the paper’s Table 2 (typical soot)
        dm_tab = [20 30 50 75 100 200 300 500 750 1000]';    % nm
        rho_tab = [1180 954 731 592 510 356 288 221 179 154]';% kg/m^3
        dp_tab = [10.1 11.7 14.0 16.1 17.8 22.7 26.1 31.3 36.0 39.8]'; % nm

        % dp -> dm via interpolation (monotonic in this range)
        dm_nm = interp1(dp_tab, dm_tab, dp_nm, 'linear', 'extrap');
        % dp -> rho via interpolation
        rho_eff = interp1(dp_tab, rho_tab, dp_nm, 'linear', 'extrap');

        % Optional: clamp outside table range to avoid wild extrapolation
        if isfield(params,'clamp') && params.clamp
            dm_nm  = max(min(dm_nm,  max(dm_tab)),  min(dm_tab));
            rho_eff= max(min(rho_eff,max(rho_tab)), min(rho_tab));
        end

    otherwise
        error('Unknown method. Use ''equations'' or ''table''.');
end
end

function v = getf(s, fname, default)
    if isfield(s, fname) && ~isempty(s.(fname)), v = s.(fname);
    else, v = default; end
end
