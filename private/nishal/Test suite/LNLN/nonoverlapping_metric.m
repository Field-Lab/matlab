function [metric,metric_su] = nonoverlapping_metric(u_spatial_log)

u_sp_norm = u_spatial_log./repelem(sqrt(sum(u_spatial_log.^2,1)),size(u_spatial_log,1),1);

%% compute pair-wise inner product

metric_su=[];
nSU = size(u_sp_norm,2);
inner_prod_mat = (u_sp_norm'*u_sp_norm);
avg_innerprod = (sum(inner_prod_mat(:)) - trace(inner_prod_mat))/(nSU^2 - nSU);

metric_su = inner_prod_mat;
metric = avg_innerprod;

end

