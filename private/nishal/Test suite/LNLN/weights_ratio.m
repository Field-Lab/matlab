function wts_log = weights_ratio(u_spatial_log,model)
u_sp_norm = u_spatial_log

%% metric for each cell ?
wts_log=[];
su_wise = zeros(model.nSU,1);
su_wise2 = zeros(model.nSU,1);
for isu=1:model.nSU
    su_vec = model.su_lowres(:,isu);
    su_vec = su_vec/norm(full(su_vec));
abcd = u_sp_norm'*su_vec;
wts = sort(abcd/max(abs(abcd)));
wts_log=[wts_log, wts(1:3)]
 
end
end