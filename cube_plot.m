function [ciao] = cube_plot (side,center)
a = -pi : pi/2 : pi;                                % Define Corners
ph = pi/4;                                          % Define Angular Orientation (‘Phase’)
x = [cos(a+ph); cos(a+ph)]*cos(ph)*side+0.5*side;
y = [sin(a+ph); sin(a+ph)]*sin(ph)*side+0.5*side;
z = [-ones(size(a)); ones(size(a))]*side/2+0.5*side;
figure
%surf(x, y, z, 'FaceColor','g')                      % Plot Cube
surf(x,y,z,'FaceAlpha',0.05)
hold on
%patch(x', y', z', 'r')                              % Make Cube Appear Solid
%hold off
axis 'equal'
%axis([ -1  1    -1  1    -1  1]*1.5)
grid on
hold on

ciao='everything ok';