function dae = AERODYN(pp, dm, rho0)
% "AERODYN" calculates the aerodynamic diameter from mass and mobility...
%   ...diameter of aggregates.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   pp: a series of 2d arrays containing primary particle data for...
%       ...different aggregates.    
%   dm: mobility diameter of aggregates
%   rho0: material density for fractal aggregates
% ----------------------------------------------------------------------- %
% 
% Outputs:
%   dae: Aerodynamic diameter of aggregates
% ----------------------------------------------------------------------- %

% constant coefficients for slip correction
A1 = 1.165;
A2 = 0.483;
A3 = 0.997;

% Cunningham correction factor (x: particle diameter)
cc = @(x) 1 + x .* (A1 + A2 * exp(-A3 ./ x));

dae0 = PAR.EQUIV(pp); % initial guess (volume-equivalent diameter)

dae = zeros(length(dae0),1); % allocate aerodynamic diameter array

m = zeros(length(dae0),1); % allocate aggregate mass array

for i = 1 : length(dae0)
    
    % calculate aggregate mass
    m(i) = rho0 * (pi/6) * sum(pp{i}(:,2) .^ 3);
    
    % define nonlinear equation
    df = @(x) x - sqrt((6/(pi*rho0)) * (m(i)/dm(i)) * (cc(dm(i))/cc(x))); 
    
    % solve the equation
    dae(i) = fzero(df, dae0(i));

end

end

