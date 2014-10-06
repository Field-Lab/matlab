function [datapoints labels means covariances] = sample_igmm_prior(n,a_0,b_0,mu_0,lambda_0,k_0,v_0)

alpha = gamrnd(a_0,b_0);
labels = sample_crp(n,alpha);

dim = length(mu_0);
num_tables = length(unique(labels));

means = zeros(num_tables,dim);
covariances = zeros(num_tables,dim,dim);
datapoints = zeros(n,dim);
dsindex = 1;
for t = 1:num_tables
    covariances(t,:,:) =iwishrnd(lambda_0,v_0);
    means(t,:) = mvnrnd(mu_0,squeeze(covariances(t,:,:))/k_0);
    num_points_at_table_t = sum(labels == t);
    datapoints(dsindex:(dsindex+num_points_at_table_t-1),:) = mvnrnd(means(t,:),squeeze(covariances(t,:,:)),num_points_at_table_t);
    dsindex = dsindex+num_points_at_table_t;
end
