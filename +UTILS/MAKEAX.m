function insetAx = MAKEAX(fig, pos, viewType)
% "MAKEAX" places a triad at a fixed position in the figure
% ----------------------------------------------------------------------- %
% Inputs:
%   fig      : handle to the figure
%   pos      : 1x4 vector [left bottom width height] in normalized units
%   viewType : '3d' or 'xy' for orientation
% ----------------------------------------------------------------------- %
% Output:
%   insetAx  : handle to the inset axes
% ----------------------------------------------------------------------- %

% create the inset axis at fixed figure coordinates
insetAx = axes('Parent', fig, 'Units', 'normalized', 'Position', pos);
hold(insetAx, 'on');
axis(insetAx, 'off');
axis(insetAx, 'equal');

% set view
switch lower(viewType)
    case '3d'
        view(insetAx, 3);
    case 'xy'
        view(insetAx, [0 90]);
    otherwise
        view(insetAx, 3);
end

% geometry
L = 2.5;
label_offset = 1.2;

% colors
xColor = hex2rgb('#BC7C7C');
yColor = hex2rgb('#B1C29E');
zColor = hex2rgb('#295F98');

% arrows
quiver3(insetAx, 0, 0, 0, L, 0, 0, 'Color', xColor, 'LineWidth', 2, ...
    'AutoScale', 'off', 'MaxHeadSize', 1);
quiver3(insetAx, 0, 0, 0, 0, L, 0, 'Color', yColor, 'LineWidth', 2, ...
    'AutoScale', 'off', 'MaxHeadSize', 1);
quiver3(insetAx, 0, 0, 0, 0, 0, L, 'Color', zColor, 'LineWidth', 2, ...
    'AutoScale', 'off', 'MaxHeadSize', 1);

% labels
xLab = [L + label_offset, 0, 0];
yLab = [0, L + label_offset, 0];
if strcmpi(viewType, 'xy')
    zLab = [0, -L * 0.6, L + 0.2];
else
    zLab = [0, 0, L + label_offset];
end

text(xLab(1), xLab(2), xLab(3), 'X', 'Color', xColor, ...
    'FontWeight', 'bold', 'FontSize', 10, 'HorizontalAlignment', 'center');
text(yLab(1), yLab(2), yLab(3), 'Y', 'Color', yColor, ...
    'FontWeight', 'bold', 'FontSize', 10, 'HorizontalAlignment', 'center');
text(zLab(1), zLab(2), zLab(3), 'Z', 'Color', zColor, ...
    'FontWeight', 'bold', 'FontSize', 10, 'HorizontalAlignment', 'center');

% visual constraints
lims = [-1 L + label_offset + 0.5];
xlim(insetAx, lims); ylim(insetAx, lims); zlim(insetAx, lims);
insetAx.DataAspectRatio = [1 1 1];
insetAx.CameraViewAngle = 6;
insetAx.CameraViewAngleMode = 'manual';

uistack(insetAx, 'top');

end
