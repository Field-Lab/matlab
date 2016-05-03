
%% smoothness estimate
list = dir('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data001/');

icell=0;
for ifile = 3:length(list)
    data = load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data001/%s',list(ifile).name),'fitGMLM_log','mask');
    nSU=4;
    filters = data.fitGMLM_log{nSU}.Linear.filter;
    mask = data.mask;
    
    u_spatial_log = zeros(sum(mask(:)),nSU);
    for isu=1:nSU
        u_spatial_log(:,isu) =filters{isu};
    end
    
    
    u_spatial_log_perm = 0*u_spatial_log;
    for idim = 1:size(u_spatial_log,1)
        perm = randperm(nSU);
        u_spatial_log_perm(idim,:)  = u_spatial_log(idim,perm);
    end
    
    [metric,metric_sus] = smoothness_metric2(u_spatial_log,mask);
    [metric_perm,metric_sus_perm] = smoothness_metric2(u_spatial_log_perm,mask);
    close all;
    [metric,metric_perm]
    icell=icell+1
    mc_data(icell).metric=metric;
    mc_data(icell).metric_perm = metric_perm;
    mc_data(icell).metric_sus = metric_sus;
    mc_data(icell).metric_sus_perm = metric_sus_perm;
end


metric=[];metric_perm=[];
for imc=1:length(mc_data)
metric = [metric;mc_data(imc).metric];
metric_perm = [metric_perm;mc_data(imc).metric_perm];
end
 
figure;
histogram(metric);
hold on;
histogram(metric_perm);
legend('Extracted sub-units','Permuted sub-units');

%% non-overlapping estimate

list = dir('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data001/');

icell=1;
for ifile = 3:length(list)
    data = load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data001/%s',list(ifile).name),'fitGMLM_log','mask');
    nSU=4;
    filters = data.fitGMLM_log{nSU}.Linear.filter;
    mask = data.mask;
    
    u_spatial_log = zeros(sum(mask(:)),nSU);
    for isu=1:nSU
        u_spatial_log(:,isu) =filters{isu};
    end
    
    
   
    [metric,metric_sus] = nonoverlapping_metric(u_spatial_log);
    
    mc_data(icell).metric=metric;
    mc_data(icell).metric_sus = metric_sus;
    mc_data(icell).u_spatial_log = u_spatial_log;
    
    u_spatial_log_perm = 0*u_spatial_log;
    for idim = 1:size(u_spatial_log,1)
        perm = randperm(nSU);
        u_spatial_log_perm(idim,:)  = u_spatial_log(idim,perm);
    end
    [metric_perm,metric_sus_perm] = nonoverlapping_metric(u_spatial_log_perm);
    
    mc_data(icell).metric_perm=metric_perm;
    mc_data(icell).metric_sus_perm = metric_sus_perm;
    mc_data(icell).u_spatial_log_perm = u_spatial_log_perm;
    
    u_spatial_log_perm2 = 0*u_spatial_log;
    for isu = 1:size(u_spatial_log,2)
        perm = randperm(size(u_spatial_log,1));
        u_spatial_log_perm2(:,isu)  = u_spatial_log(perm,isu);
    end
    [metric_perm2,metric_sus_perm2] = nonoverlapping_metric(u_spatial_log_perm2);
    
    mc_data(icell).metric_perm2=metric_perm2;
    mc_data(icell).metric_sus_perm2 = metric_sus_perm2;
    mc_data(icell).u_spatial_log_perm2 = u_spatial_log_perm2;

    
    u_spatial_log_randn = randn(size(u_spatial_log,1),size(u_spatial_log,2));
    [metric_randn,metric_sus_randn] = nonoverlapping_metric(u_spatial_log_randn);
    
    mc_data(icell).metric_randn=metric_randn;
    mc_data(icell).metric_sus_randn = metric_sus_randn;
    mc_data(icell).u_spatial_log_randn = u_spatial_log_randn;

    u_spatial_log_rand = rand(size(u_spatial_log,1),size(u_spatial_log,2));
    [metric_rand,metric_sus_rand] = nonoverlapping_metric(u_spatial_log_rand);
    
    mc_data(icell).metric_rand=metric_rand;
    mc_data(icell).metric_sus_rand = metric_sus_rand;
    mc_data(icell).u_spatial_log_rand = u_spatial_log_rand;
    mc_data(icell).filename = list(ifile).name;
    
    close all;
    [metric,metric_perm,metric_perm2,metric_randn,metric_rand]
    icell=icell+1
    
end


metric=[];metric_perm=[];metric_perm2=[];metric_randn=[];metric_rand=[];
for imc=1:length(mc_data)
metric = [metric;mc_data(imc).metric];
metric_perm = [metric_perm;mc_data(imc).metric_perm];
metric_perm2 = [metric_perm2;mc_data(imc).metric_perm2];
metric_randn = [metric_randn;mc_data(imc).metric_randn];
metric_rand = [metric_rand;mc_data(imc).metric_rand];
end
 
figure;
histogram(metric);
hold on;
histogram(metric_perm);
hold on;
histogram(metric_perm2,20);
hold on;
histogram(metric_randn);
hold on;
histogram(metric_rand);
legend('Extracted sub-units','Permuted sub-units su wise','Permuted sub-units dim wise','randn','rand');

% sparsity of sub-units 
u_l=[];
for imc=1:length(mc_data)
u_l = [u_l;mc_data(imc).u_spatial_log(:)];
end

u_l=[];
for imc=1:length(mc_data)
    u_spatial_log = mc_data(imc).u_spatial_log;
u_sp_norm = u_spatial_log./repelem(sqrt(sum(u_spatial_log.^2,1)),size(u_spatial_log,1),1);
u_l = [u_l;u_sp_norm(:)];
end
figure;
histogram(u_l(:));

%% non-overlapping index 2

list = dir('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data001/');

for imc = 1:length(mc_data)
    imc
    data = load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data001/%s',mc_data(imc).filename),'mask');

    [self,cross] = nonoverlapping_metric2(mc_data(imc).u_spatial_log,data.mask);
    mc_data(imc).self_xcorr = self;
    mc_data(imc).cross_corr = cross;
     
end

sum_xcorr = zeros(159,79);
sum_self_corr = sum_xcorr;

for imc=1:length(mc_data)
sum_xcorr = mc_data(imc).cross_corr;
sum_self_corr = mc_data(imc).self_xcorr;
end

figure;
contourf(sum_xcorr(75:85,35:45));
figure;
contourf(sum_self_corr(75:85,35:45));

%% non-overlapping index 3

list = dir('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data001/');

for imc = 1:length(mc_data)
    imc
    data = load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data001/%s',mc_data(imc).filename),'mask');

    [fits,clog] = nonoverlapping_metric3(mc_data(imc).u_spatial_log,data.mask);
    mc_data(imc).fits = fits;
    mc_data(imc).clog = clog;
     
end

% all pairs

cl=[];
for imc=1:length(mc_data)
cl = [cl;mc_data(imc).clog.cc_ratio];    
end
cl_d = -2*log( cl);
cl_dist = sqrt(cl_d(cl_d>0));
figure;
histogram(real(cl_dist));
set(gca,'yTick',[]);

% nearest neighbors
figure
for nn=1:3
cl = [];

for imc=1:length(mc_data)
    for isu=1:length(mc_data(imc).clog.center_log2)
        center = mc_data(imc).clog.center_log2(isu);
        sort_centers =  sort(mc_data(imc).clog.other_centers(isu,:),'descend');
        nearest_center = sort_centers(nn);
        cl = [cl;nearest_center/center];
    end
end

cl_d = -2*log( cl);
cl_dist = sqrt(cl_d(cl_d>0));
hold on;
histogram(real(cl_dist));
end

% plot 2D gaussians
[X,Y] = meshgrid([-4:0.1:4],[-4:0.1:5]);
figure;
Z = normpdf(X,0,1).*normpdf(Y,0,1);
surf(X,Y,Z);
hold on;
Z = normpdf(X,0,1).*normpdf(Y,1.5,1);
surf(X,Y,Z);
set(gca,'visible','off')