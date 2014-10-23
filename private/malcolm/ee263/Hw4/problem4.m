m = 20;
D = 12;
n = length(h);

% make convoluton matrix H
H = nan(n+m-1,m);
for i = 1:n+m-1
    head = max(i-n,0);
    tail = max(m-i,0);
    H(i,:) = [zeros(1,head) h(min(i,n):-1:1-min(m-n-head,0))' zeros(1,tail)];
end

HminusD = H([1:D D+2:n+m-1],:);
HD = H(D+1,:);

A = [HminusD'*HminusD HD'; HD 0];
glambda = A\[zeros(m,1) ; 1];
g = glambda(1:m);

% make plots
subplot(3,1,1)
plot(0:n-1,h)
subplot(3,1,2)
plot(0:m-1,g)
subplot(3,1,3)
plot(0:n+m-2,H*g)

% compute z by convolving g and y
z = conv(g, y);
t = 1:length(y);
z = z(t);
figure; 
subplot(3,1,1)
plot(t,y,t,z,'--')
legend({'y','z'})
subplot(3,1,2)
hist(y)
subplot(3,1,3)
hist(z(D+1:end-D))