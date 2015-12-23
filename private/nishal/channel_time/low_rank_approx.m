
data_cut=data(:,1:600);
cvx_begin
variable X(512,600)

norm_nuc(X)<=100

 minimize (norm(X(:)-data_cut(:)))
cvx_end

%% 
mean_dat=mean(data,2);
data_cut=data(:,1:1000)-repmat(mean_dat,[1,1000]);
[U,S,V]=svd(data_cut);
