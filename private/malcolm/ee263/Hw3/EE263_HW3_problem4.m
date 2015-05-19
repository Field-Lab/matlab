t = 1:N;
A = [t' repmat(eye(24),7,1)];
ap = A\y;
yhat = A*ap;
subplot(2,1,1);
plot(t,y,'o',t,yhat,'x');

t2 = 169:192;
ytomhat = [t2' eye(24)]*ap;
rms = sqrt(sum((ytomhat-ytom).^2)/24);
subplot(2,1,2);
plot(t2, ytom, 'o', t2, ytomhat, 'x');