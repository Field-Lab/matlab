function drawEllipse(mu,V,color)
% Draws a simple ellipse

theta = 2*pi*[0:0.01:1];
x = [cos(theta); sin(theta)];
[V,D] = eig(inv(V));
x = V*sqrt(D)*x;

for ii=1:length(mu)
    x(ii,:) = x(ii,:) + mu(ii);
end

if nargin<3
    plot(x(1,:),x(2,:));
else
    plot(x(1,:),x(2,:),color);
end

end % drawEllipse
