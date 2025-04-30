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
    'red',     [1 0 0], ...
    'green',   [0 1 0], ...
    'blue',    [189, 221, 228]/255, ...
    'purple',  [205, 193, 255]/255, ...
    'yellow',  [1 1 0], ...
    'orange',  [243, 158, 96]/255, ...
    'pink',    [1 0.4 0.7], ...
    'cyan',    [0 1 1], ...
    'gray',    [0.5 0.5 0.5]);

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