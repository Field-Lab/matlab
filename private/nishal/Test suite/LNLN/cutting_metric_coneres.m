function [metric_comb,su_wise_comb] = cutting_metric_coneres(model,u_spatial_log)

u_sp_norm = u_spatial_log./repelem(sqrt(sum(u_spatial_log.^2,1)),size(u_spatial_log,1),1);

%% metric for each cell ?
su_wise = zeros(model.nSU,1);
su_wise2 = zeros(model.nSU,1);
for isu=1:model.nSU
    su_vec = double(model.cone_su_idx==isu);
    su_vec = su_vec/norm(full(su_vec));
abcd = u_sp_norm'*su_vec;
%su_wise(isu) = sqrt(max(abcd.^2))/norm(abcd);
su_wise(isu) = sqrt((sum(abcd.^2) - max(abcd).^2))/norm(abcd);
su_wise2(isu) = max(abs(abcd));
end

su_wise_comb = [su_wise,su_wise2];

%% overall metric
metric = mean(su_wise);
metric2 =mean(su_wise2);
metric_comb = [metric,metric2];

end