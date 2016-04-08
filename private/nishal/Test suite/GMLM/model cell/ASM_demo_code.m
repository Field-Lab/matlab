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
nsamples = 100000;
nsample_batch=500;
x = randn(2,nsamples);
k1 = [1,0.2];
k2 = [0.2,1];

%f = @(x) (x.*(x>0)).^2 + 0.001);
%fd= @(x) 2*(x.*(x>0));

f = @(x) exp(x);
fd = @(x) exp(x);
lam = (f(k1*x) + f(k2*x))/10;
spk =poissrnd(lam);

k1_0 = [0.8,0.2];k2_0=[0.3,1];

xcnt=0;
ll=[];ll_ub_full=[];ll_ub_batch=[];
range = -2:0.1:2; k1_fixed=0.2;k2_fixed=1;
for k1_ctr=range
    xcnt=xcnt+1;
    ycnt=0;
    for k2_ctr=range
        ycnt=ycnt+1;
        
        k1_grid = [k1_ctr,k1_fixed];
        k2_grid = [k2_ctr,k2_fixed];
        
        lam_grid =  (f(k1_grid*x) + f(k2_grid*x))/10;
        ll(xcnt,ycnt) = (sum(lam_grid) - spk*log(lam_grid'))/nsamples;
        
        lam_approx = (f(k1_0*x) + f(k2_0*x))/10;
        ll_ub_full(xcnt,ycnt) = (sum(lam_grid) - (spk*log(lam_approx') + ... 
        spk*((1./lam_approx).* (fd(k1_0*x).*((k1_grid-k1_0)*x/10) + fd(k2_0*x).*((k2_grid-k2_0)*x/10) )  )' ))/nsamples;
        
        x_batch = x(:,1:nsample_batch);spk_batch = spk(1:nsample_batch);
        lam_approx_batch = (f(k1_0*x_batch) + f(k2_0*x_batch))/10;
        lam_grid_batch =  (f(k1_grid*x_batch) + f(k2_grid*x_batch))/10;
        ll_ub_batch(xcnt,ycnt) = (sum(lam_grid_batch) - (spk_batch*log(lam_approx_batch') + ... 
        spk_batch*((1./lam_approx_batch).* (fd(k1_0*x_batch).*((k1_grid-k1_0)*x_batch/10) + fd(k2_0*x_batch).*((k2_grid-k2_0)*x_batch/10) )  )' ))/nsample_batch;
        
    end
end
npts = length(range);
figure;
contourf(repelem(range',1,npts),repelem(range,npts,1),ll);


figure;
contourf(repelem(range',1,npts),repelem(range,npts,1),ll_ub_full);

figure;
contourf(repelem(range',1,npts),repelem(range,npts,1),abs(ll_ub_full-ll));
hold on;
plot(k1_0(1),k2_0(1),'r.','MarkerSize',20);
caxis([-1.5,1.5]);

figure;
contourf(repelem(range',1,npts),repelem(range,npts,1),(ll_ub_batch-ll));
hold on;
plot(k1_0(1),k2_0(1),'r.','MarkerSize',20);
caxis([-1.5,1.5]);


figure;
contourf(repelem(range',1,npts),repelem(range,npts,1),(ll_ub_batch-ll_ub_full),100);
hold on;
plot(k1_0(1),k2_0(1),'r.','MarkerSize',20);


%% plot distance metric
nsamples = 100000;
nsample_batch=500;
x = randn(2,nsamples);
x_batch = x(:,1:nsample_batch);spk_batch = spk(1:nsample_batch);

k1 = [1,0.2];
k2 = [0.2,1];

f = @(x) exp(x);
fd = @(x) exp(x);
lam = (f(k1*x) + f(k2*x))/10;
spk =poissrnd(lam);

k1_0 = [0.8,0.2];k2_0=[0.3,1];

k1_ref = [1,0.2];k2_ref = [0.2,1];
xcnt=0;
dist=[];
range = -2:0.1:2; k1_fixed=0.2;k2_fixed=1;
x_spk = x(:,spk>0);
lam_ref =  (f(k1_ref*x_spk) + f(k2_ref*x_spk))/10;
weight_ref = [fd(k1_ref*x_spk)/10;fd(k2_ref*x_spk)/10]./repelem(lam_ref,2,1);
        
for k1_ctr=range
    xcnt=xcnt+1;
    ycnt=0;
    for k2_ctr=range
        ycnt=ycnt+1;
        
        k1_grid = [k1_ctr,k1_fixed];
        k2_grid = [k2_ctr,k2_fixed];
        
        lam_grid =  (f(k1_grid*x_spk) + f(k2_grid*x_spk))/10;
        weight_grid = [fd(k1_grid*x_spk)/10;fd(k2_grid*x_spk)/10]./repelem(lam_grid,2,1);
         
        dist(xcnt,ycnt) = mean(sum(abs(weight_ref-weight_grid),1),2);
    end
end

figure;
contourf(repelem(range',1,npts),repelem(range,npts,1),dist);
hold on;
plot(k1_ref(1),k2_ref(1),'r.','MarkerSize',20);
%caxis([-1.5,1.5]);

