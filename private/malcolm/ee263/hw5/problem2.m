% problem 2b

A = [0 1; -1 0];
[x,y] = meshgrid(-2:2,-2:2);
figure();
quiver(x,y,y,-x);
title('Problem 2b');

% problem 2c
% plot some circular orbits
x0 = [0.5 1 2;0.5 1 2];
t = 0:0.1:2*pi;
x = nan(2,length(t),size(x0,2));
for i = 1:length(t)
    for j = 1:size(x0,2);
        x(:,i,j) = expm(t(i)*A)*x0(:,j);
    end
end
figure();
plot(x(1,:,1),x(2,:,1),'.',x(1,:,2),x(2,:,2),'.',x(1,:,3),x(2,:,3),'.');
title('Problem 2c');