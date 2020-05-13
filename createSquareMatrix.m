function [ cartesian_matrix, theta_matrix, phi_matrix ] = createSquareMatrix( c1, c2, ~, c4, nm, plot3D)
% Creates square computing mesh in cartesian coordinate system, and
% determines computing mesh points in cartesian and spherical coordinate
% system covering the object area. Computing mesh can be
% utilized e.g. with geographical or optical simulations
% 
% c1, c2 .. define four (3) corners in 3D space [x y z] in rotational order
% nm:       total number of computing points [n m]
% plot3D:   1 to plot computing mesh and detector area. 
% 
% Antti Järvinen 11.5.2020

% An example:
% [ cartesian_matrix, theta_matrix, phi_matrix ] = createSquareMatrix( [2 2 -2], [2 2 2], [2 -2 2], [2 -2 -2], [10 10], 1)


% Discrete detector matrix points (n number) in cartesian coordinate system
x = zeros(1,nm(1)*nm(2));
y = zeros(1,nm(1)*nm(2));
z = zeros(1,nm(1)*nm(2));

% Use c1 as a starting corner. Order of corners is clockwise, when looking
% back side of the detector (looking through the detector to direction of
% sphere)

% Define c1 -> c2 / (n-1), and c1 -> c4 / (m-1) [xi, yi, zi] vectors
c2c1 = c2 - c1;
c2c1n = c2c1 / (nm(1)-1);
c4c1 = c4 - c1;
c4c1n = c4c1 / (nm(2)-1);

% Calculate each point [x y x] of detector matrix in 3D space;
counter = 1;
cartesian_matrix = zeros(nm(1),nm(2),3); 
for i = 0:(nm(1)-1)
    for j = 0:(nm(2)-1)
        
        % Calculate absolute [x y x] positions
        position = c1 + (c2c1n * i) + (c4c1n * j);
        x(counter) = position(1);
        y(counter) = position(2);
        z(counter) = position(3);
        
        % Unit vector = [x y z] / norm([x y x])
        euclidean_norm = norm([x(counter) y(counter) z(counter)]);
        xyz_unit = [x(counter) y(counter) z(counter)] / euclidean_norm;
        
        % Fill detector matrix with unit vector values
        cartesian_matrix(i+1,j+1,1) = xyz_unit(1);
        cartesian_matrix(i+1,j+1,2) = xyz_unit(2);
        cartesian_matrix(i+1,j+1,3) = xyz_unit(3);
        
        counter = counter + 1;
    end
end

% Convert cartesian [x y z] coordinates to spherical coordinates (theta,
% phi, r)
phi_matrix = zeros(nm(1),nm(2));
theta_matrix = zeros(nm(1),nm(2));
for i = 1:nm(1)
    for j = 1:nm(2)
        % azimuth(phi) = atan2(y,x)
        phi_matrix(i,j) = atan2(detector_matrix_car(i,j,2),detector_matrix_car(i,j,1));
        % elevation(theta) = atan2(z,sqrt(x.^2 + y.^2))
        theta_matrix(i,j) = atan2(detector_matrix_car(i,j,3),sqrt(detector_matrix_car(i,j,1).^2 + detector_matrix_car(i,j,2).^2));
        % r = sqrt(x.^2 + y.^2 + z.^2)
        % Euclidean unit vectors (r=1)
    end
end

%% Plots computing mesh points in 3D space. Points located in the detector area are drawn in different color.
if plot3D == 1 

    % Plots the computing mesh. 
    figure(10) 
    
    % Plot absolute position of detector in 3D space
    scatter3(x,y,z, 'Marker', '.','MarkerEdgeColor','k','MarkerFaceColor','k') 
    hold on 
    
    % Plot detector as unit vectors
    xu = cartesian_matrix(:,:,1);
    yu = cartesian_matrix(:,:,2);
    zu = cartesian_matrix(:,:,3);
    scatter3((xu(:)'),(yu(:)'),(zu(:)'), 'Marker', '.','MarkerEdgeColor','r','MarkerFaceColor','r') 
    
    % Set labels
    xlabel('x') 
    ylabel('y') 
    zlabel('z') 

    % Plots axis through origin. 
    plot3([0 0],[0 0],[-1 1],'k','linewidth',2) 
    plot3([0 0],[-1 1],[0 0],'k','linewidth',2) 
    plot3([-1 1],[0 0],[0 0],'k','linewidth',2) 
    
    % Axis limits 
    xlim([min([-2 min(x)]) max([2 max(x)])]) 
    ylim([min([-2 min(y)]) max([2 max(y)])])  
    zlim([min([-2 min(z)]) max([2 max(z)])]) 
end
