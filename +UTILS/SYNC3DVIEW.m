function SYNC3DVIEW(axSource, axTarget)
% "SYNC3DVIEW" Synchronize 3D view and zoom between two axes
%   sync3DView(axSource, axTarget) copies the view orientation, axis
%   limits, and full camera settings from one 3D axes (axSource) to another
%   (axTarget), ensuring both display the same orientation and spatial
%   extent. This is useful for tiledlayout figures with multiple 3D plots.
% ----------------------------------------------------------------------- %
%   INPUTS:
%       axSource - Handle to the source axes object (must contain a 3D plot)
%       axTarget - Handle to the target axes object to be updated
% ----------------------------------------------------------------------- %
%   OUTPUTS:
%       None. The function modifies axTarget in place to match axSource.
% ----------------------------------------------------------------------- %
%   EXAMPLE:
%       t = tiledlayout(1,2);
%       nexttile; ax1 = gca; surf(peaks);
%       nexttile; ax2 = gca; surf(peaks);
%       sync3DView(ax1, ax2);  % Make both views identical
% ----------------------------------------------------------------------- %

   % Safely copy azimuth/elevation view vector
    azel = get(axSource, 'View'); % returns a 1x2 vector
    view(axTarget, azel);         % apply to target

    % Copy axis limits
    xlim(axTarget, xlim(axSource));
    ylim(axTarget, ylim(axSource));
    zlim(axTarget, zlim(axSource));

    % % Copy full camera properties
    % campos(axTarget, campos(axSource));
    % camtarget(axTarget, camtarget(axSource));
    % camup(axTarget, camup(axSource));
    % camva(axTarget, camva(axSource));
    
end
