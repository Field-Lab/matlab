
% make on and off parasol mosaic

on = model_population_stas('coneLatticeOrientation',0/3,'gridsz',256);
off = model_population_stas('coneLatticeOrientation',pi/6,'gridsz',256);

% plot cells
n=100;
angle = 0:2*pi/n:2*pi;            % vector of angles at which points are drawn
R = 3*on.coneGaussSd;                         % Unit radius
x = R*cos(angle);  y = R*sin(angle);   % Coordinates of the circle

figure;
plot(on.conesX,on.conesY,'r*');
hold on;

for icone=1:on.nCones
    plot(x+on.conesX(icone),y+on.conesY(icone),'r');                      % Plot the circle 
    axis equal;
    grid on;
end
hold on;

R = 3*off.coneGaussSd;                         % Unit radius
x = R*cos(angle);  y = R*sin(angle);   % Coordinates of the circle

plot(off.conesX,off.conesY,'b*');
hold on;
n=100;
for icone=1:off.nCones
    plot(x+off.conesX(icone),y+off.conesY(icone),'b');                      % Plot the circle 
    axis equal;
    grid on;
end
