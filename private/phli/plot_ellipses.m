function plot_ellipses(cell_array_of_structs, color)

if nargin < 2
    color = 'b';
end

circle_samples = 0:0.05:2*pi;
x_circle = cos(circle_samples);
y_circle = sin(circle_samples);

%figure
hold on;

for i=1:length(cell_array_of_structs)
    % grab parameters
    sd = cell_array_of_structs{i}.sd;
    mn = cell_array_of_structs{i}.mean;
    angle = cell_array_of_structs{i}.angle;
    
    % rotate by angle and stretch
    angle = -angle;
    R = rmatrix2(angle / (2*pi) * 360);
    L = [sd(1) 0; 0 sd(2)];
    
    % Set coordinates and close loop
    X = R * L * [x_circle; y_circle];
    X(:,end+1) = X(:,1);
    
    % Offset to center points
    X(1,:)=(X(1,:)+mn(1));
    X(2,:)=(abs((X(2,:)+mn(2))));
    
    plot(X(1,:), X(2,:), color);
end
