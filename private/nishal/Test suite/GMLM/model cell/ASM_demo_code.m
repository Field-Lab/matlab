x = randn(2,10000);
k1 = [1,0.2];
k2 = [0.2,1];

%f = @(x) (x.*(x>0)).^2 + 0.001);
%fd= @(x) 2*(x.*(x>0));

f = @(x) exp(x);
fd = @(x) exp(x);
lam = (f(k1*x) + f(k2*x))/10;
spk =poissrnd(lam);

x_spk = x(:,spk>0);

w1 = fd(k1*x)./(f(k1*x) + f(k2*x));
w2 = fd(k2*x)./(f(k1*x) + f(k2*x));

figure;
plot(x(1,:),x(2,:),'k.');
hold on;
plot([-5,5],[0,0],'g');
hold on;
plot([0,0],[-5,5],'g');
hold on;
scatter(x(1,spk>0),x(2,spk>0),15,w1(spk>0),'filled');colormap hot
sta = (x*spk'/sum(spk));
hold on;
plot([0,3*sta(1)],[0,3*sta(2)],'LineWidth',2);
hold on;
plot([0,5*k1(1)],[0,5*k1(2)],'b','LineWidth',2);
hold on;
plot([0,5*k2(1)],[0,5*k2(2)],'b','LineWidth',2);


%% Algorithm 
k_init1 = randn(1,2);k_init2=randn(1,2);
k_iter1=k_init1;k_iter2=k_init2;

for iter=1:20
w1 = fd(k_iter1*x)./(f(k_iter1*x) + f(k_iter2*x));
w2 = fd(k_iter2*x)./(f(k_iter1*x) + f(k_iter2*x));

figure;
plot(x(1,:),x(2,:),'k.');
hold on;
scatter(x(1,spk>0),x(2,spk>0),15,w1(spk>0),'filled');colormap hot
hold on;
plot([-1,1],[0,0],'g');
hold on;
plot([0,0],[-1,1],'g');

sta = (x*spk'/sum(spk));
hold on;
plot([0,3*sta(1)],[0,3*sta(2)],'LineWidth',2);
hold on;
plot([0,5*k1(1)],[0,5*k1(2)],'b','LineWidth',2);
hold on;
plot([0,5*k2(1)],[0,5*k2(2)],'b','LineWidth',2);

hold on;
plot([0,3*k_iter1(1)],[0,3*k_iter1(2)],'m','LineWidth',2);
hold on;
plot([0,3*k_iter2(1)],[0,3*k_iter2(2)],'m','LineWidth',2);

k_iter1 = (x*(w1.*spk)'/sum(w1.*spk));k_iter1=k_iter1';
k_iter2 = (x*(w2.*spk)'/sum(w2.*spk));k_iter2=k_iter2';
title(sprintf('iter: %d',iter));
end

%% True log-likelihood surface - a LOT of data!
nsamples = 10000;
x = randn(2,10000);
k1 = [1,0.2];
k2 = [0.2,1];

%f = @(x) (x.*(x>0)).^2 + 0.001);
%fd= @(x) 2*(x.*(x>0));

f = @(x) exp(x);
fd = @(x) exp(x);
lam = (f(k1*x) + f(k2*x))/10;
spk =poissrnd(lam);

xcnt=0;
ll=[];
range = -1:0.1:1;
for k1_grid=-1:0.1:1
    xcnt=xcnt+1;
    ycnt=0;
    for k2_grid=-1:0.1:1
        ycnt=ycnt+1;
        lam_grid =  sum((f(k1_grid*x) + f(k2_grid*x))/10);
        ll(xcnt,ycnt) = (sum(lam_grid) - spk*log(lam_grid'))/nsamples;
        
    end
end
npts = length(range);
figure;
contourf(repelem(range',1,npts),repelem(range,npts,1),ll);

% get approximation at a particular point! k0!! 


