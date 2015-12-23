nmax = 10;
ebar = zeros(nmax,1);
ebarcv = zeros(nmax,1);

for n = 1:nmax
    G = zeros(N-n,n+1);
    for t = (n+1):N
        G(t-n,:) = [u(t) y((t-1):-1:(t-n))'];
    end
    ysub = y((n+1):N);
    ab = G\ysub;
    yhat = G*ab;
    ebar(n) = sum((yhat-ysub).^2)/sum(ysub.^2);
    
    Gcv = zeros(N-n,n+1);
    for t = (n+1):N
        Gcv(t-n,:) = [ucv(t) ycv((t-1):-1:(t-n))'];
    end
    ysubcv = ycv((n+1):N);
    yhatcv = Gcv*ab;
    ebarcv(n) = sum((yhatcv-ysubcv).^2)/sum(ysubcv.^2);
end

subplot(2,1,1)
plot(1:nmax, ebar)
subplot(2,1,2)
plot(1:nmax, ebarcv)