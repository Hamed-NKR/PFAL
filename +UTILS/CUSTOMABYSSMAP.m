function cmap = CUSTOMABYSSMAP(targetColor, n)
% "CUSTOMABYSSMAP" Generates a custom colormap with abyss-like gradients,
%   replacing blue with the target color.
% ----------------------------------------------------------------------- %
%   Inputs:
%       targetColor - String: 'red', 'green', 'purple', 'yellow', 'orange', etc.
%       n           - Number of colormap levels (default = 256)
% ----------------------------------------------------------------------- %
%   Output:
%       cmap        - n-by-3 colormap array
% ----------------------------------------------------------------------- %

if nargin < 2
    n = 256;
end

% Load original 'abyss' colormap (requires cmocean)
orig_map = abyss(n);

% Define RGB triplets for common target colors
color_dict = struct( ...
    'red',     hex2rgb('#DA6C6C'),...
    'green',   hex2rgb('#659287'),...
    'blue',    hex2rgb('#6096B4'),...
    'purple',  hex2rgb('#756AB6'),...
    'yellow',  hex2rgb('#F39E60'),...
    'orange',  hex2rgb('#E16A54'),...
    'pink',    hex2rgb('#E5989B'),...
    'cyan',    hex2rgb('#9ECAD6'),...
    'gray',    hex2rgb('#98A1BC'));

% Normalize target color name
targetColor = lower(targetColor);

% Check if color exists
if ~isfield(color_dict, targetColor)
    error('Unknown target color "%s". Supported colors: %s', ...
        targetColor, strjoin(fieldnames(color_dict), ', '));
end

target_rgb = color_dict.(targetColor);

% Original blue channel to be replaced
blue_gradient = orig_map(:,3);  % assuming blue dominates in 'abyss'

% Create new colormap
cmap = orig_map;
cmap(:,1) = target_rgb(1) * blue_gradient;
cmap(:,2) = target_rgb(2) * blue_gradient;
cmap(:,3) = target_rgb(3) * blue_gradient;

end