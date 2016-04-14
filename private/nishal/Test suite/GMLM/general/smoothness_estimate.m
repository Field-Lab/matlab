
list = dir('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data001/');

icell=0;
for ifile = 81:length(list)
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
