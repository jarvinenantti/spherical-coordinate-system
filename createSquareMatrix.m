function [ detector_matrix ] = createSquareMatrix( c1, c2, ~, c4, nm, plot3D)
% Creates square mesh, and transforms it to computing mesh in a surface of
% sphere having Euclidean unit distance to origo. Computing mesh can be
% utilized e.g. with geographical or optical simulations
% 
% c1, c2 .. define four (3) corners in 3D space [x y z] in rotational order
% nm:       total number of computing points [n m]
% plot3D:   1 to plot computing mesh and detector area. 
% 
% Antti Järvinen 11.5.2020

% An example:
% [ detector_matrix ] = createSquareMatrix( [2 2 -2], [2 2 2], [2 -2 2], [2 -2 -2], [10 10], 1)


% Discrete detector matrix points (n number) in 3D space
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
detector_matrix = zeros(nm(1),nm(2),3); 
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
        detector_matrix(i+1,j+1,1) = xyz_unit(1);
        detector_matrix(i+1,j+1,2) = xyz_unit(2);
        detector_matrix(i+1,j+1,3) = xyz_unit(3);
        
        counter = counter + 1;
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
    xu = detector_matrix(:,:,1);
    yu = detector_matrix(:,:,2);
    zu = detector_matrix(:,:,3);
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
    xlim([-2.5 2.5]) 
    ylim([-2.5 2.5]) 
    zlim([-2.5 2.5]) 
end
