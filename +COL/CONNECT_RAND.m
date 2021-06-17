function pp2 = CONNECT_RAND(pp1, pp2)
% "CONNECT_RAND" sticks two approaching aggregates to each other in a...
    % ...random manner.
% ----------------------------------------------------------------------- %

% Inputs/Outputs:
    % pp1: Primary particle info array for aggregate 1
    % pp2: ~ for aggregate 2 (This is the one that moves with respect to...
        % ...the other)
% ----------------------------------------------------------------------- %

n_pp1 = size(pp1, 1); % Number of primary particles within aggregate 1
n_pp2 = size(pp2, 1); % ~ aggregate 2
pp1_temp = repelem(pp1, n_pp2, 1);
pp2_temp = repmat(pp2, n_pp1, 1);

dist_pp = sqrt(sum((pp1_temp(:,3:5) - pp2_temp(:,3:5)).^2, 2));
    % List of distances of primary particle centers
ind_cnt = find(dist_pp == min(dist_pp), 1); % Finding pairs with minimum...
    % ...distance
ind_cnt = [ind_cnt, pp1_temp(ind_cnt,1), pp2_temp(ind_cnt,1)]; % Adding...
    % ...the index numbers of the nearest pair
ind_cnt = [ind_cnt, find(pp1(:,1) == ind_cnt(2)),...
    find(pp2(:,1) == ind_cnt(3))];

chk = 1; % Overlap checking criterion   
while ~ isempty(find(chk == 1, 1))
    
    pp2_com = COL.EQUIV({pp2}, n_pp2); % Aggregate 2's center of mass
    pp2_com = repmat(pp2_com, n_pp2, 1);
    
    % Pre-deposition intrinsic rotation of aggregate
    angs = 2 * pi * rand(3,1); % A set of 3 Euler angles (yaw, pitch,...
        % ...and roll)
    yaw = [cos(angs(1)), -sin(angs(1)), 0;...
        sin(angs(1)), cos(angs(1)), 0; 0, 0, 1];
        % transformation matrix for yaw rotation
    pitch = [cos(angs(2)), 0, sin(angs(2));...
         0, 1, 0; -sin(angs(2)), 0, cos(angs(2))];
        % transformation matrix for pitch rotation
    roll = [1, 0, 0; 0, cos(angs(3)), -sin(angs(3));...
        0, sin(angs(3)), cos(angs(3))];
        % transformation matrix for roll rotation
    rot = yaw * pitch * roll; % The net rotation matrix
    % Obtaining the post rotation coordinates of primaries
    pp2(:,3:5) = (rot * (pp2(:,3:5) - pp2_com)')' + pp2_com;
    
    phi = 2 * pi * rand; % A random azimuthal orientation for deposition...
        % ...of aggregate 2 on 1
    theta = pi * rand; % A random polar orientation for deposition

    % Spherical coordinates for the aggregate 2's contact primary...
        % ...particle center
    disp1 = [abs(pp1(ind_cnt(4), 2) + pp2(ind_cnt(5), 2)) / 2;...
        theta; phi];
    
    % Calculating the new center position for the aggregate 2's contact...
        % ...primary particle
    dr_pp2 = pp1(ind_cnt(4), 3:5) + ...
        [disp1(1) * sin(disp1(2)) * cos(disp1(3)), ...
        disp1(1) * sin(disp1(2)) * sin(disp1(3)), ...
        disp1(1) * cos(disp1(2))] - pp2(ind_cnt(5), 3:5);
    pp2(:, 3:5) = pp2(:, 3:5) + repmat(dr_pp2, n_pp2, 1);
    
    % Checking for overlapping of two aggregates
    pp2_temp = repmat(pp2, n_pp1, 1);
    chk = COL.OVR([pp1_temp(:, 3:5), pp2_temp(:, 3:5)],...
        [pp1_temp(:,2), pp2_temp(:,2)]);
    
end

end

