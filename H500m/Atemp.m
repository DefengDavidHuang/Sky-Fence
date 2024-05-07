% Generate sample data
t = linspace(0, 10, 100); % Time vector
y = linspace(-1, 1, 10); % Set of points
[x, y] = meshgrid(t, y); % Create a grid of time vs. points
z = sin(x) .* cos(y); % Function of time and points

% Create 3D plot
figure;
plot3(x, y, z, 'b', 'LineWidth', 2);
xlabel('Time');
ylabel('Points');
zlabel('Z');
title('3D Plot: Z = sin(Time) * cos(Points)');
grid on;