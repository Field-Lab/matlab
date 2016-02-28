function plot_cells_elecs(on,elecs,icell)

% plot cells
n=100;
angle = 0:2*pi/n:2*pi;            % vector of angles at which points are drawn
R = 2*on.coneGaussSd;                         % Unit radius
x = R*cos(angle);  y = R*sin(angle);   % Coordinates of the circle


%figure;
%plot(on.conesX,on.conesY,'r*');
hold on;

for icone=1:on.nCones
    plot(x+on.conesX(icone),y+on.conesY(icone),'k');                      % Plot the circle 
    axis equal;
    grid on;
end
hold on;

hold on;
plot(elecs.x,elecs.y,'b.')

%hold on;
%icell=15;

hold on;
plot(elecs.x(on.elecs.weight_elecs(icell,:)>0),elecs.y(on.elecs.weight_elecs(icell,:)>0),'r.');

end