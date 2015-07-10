function [fitGLM,f_val] = fitGMLM_EM_bias(binnedResponses,mov,filteredStimDim2,nFrontEnds2,interval)

%% Works for non-gaussian stimuli and there is bias inside exponential 
 
%% Data
global_vars_GMLM_afterSTC
 binnedResponses_global=binnedResponses; 
 filteredStimDim=filteredStimDim2;
 nFrontEnds=nFrontEnds2;
 
mov_filtered=mov;
phase2=0;

%% Initialization.
% 
%  Random initialization
initialFilters = (2/nFrontEnds)*(rand(filteredStimDim*nFrontEnds,1)-0.5);
initialFilters = 2*(rand(filteredStimDim*nFrontEnds,1)-0.5);
 % Add bias in initial filter
% initialFilters=[initialFilters;0];
 
% 
% % k-means of spike triggered stimuli initialization
% tsp = find(binnedResponses_global~=0);
% mov_spiked = mov_filtered(:,tsp);
% [label,centroids] = kmeans(mov_spiked',nFrontEnds);
% centroids=centroids';
% initialFilters = centroids(:);

% Extract filters 
x=initialFilters;

filters=zeros(nFrontEnds,filteredStimDim);
bias=zeros(nFrontEnds,1);
for ifilter=1:nFrontEnds
filters(ifilter,:)=x((ifilter-1)*filteredStimDim+1:ifilter*filteredStimDim);
bias(ifilter)=2*(rand(1,1)-0.5);
end
%mu=sum(binnedResponses)/(((nFrontEnds +1)^2)* size(mov_filtered,2)*interval/120);
mu=0;
%% E-M MEL . 
N =size(mov_filtered,2);
dt = interval/120;
% calculate alphas 
% SU outputs
togo=1;
f_val=Inf;
tol=1e-5;
% make 0 mean. 
mov_filtered_b =gpuArray([mov_filtered;ones(1,size(mov_filtered,2))]);
f_val_log=gpuArray([]); icnt=0;

mov_filtered = gpuArray(mov_filtered);
binnedResponses_global = gpuArray(binnedResponses_global);
tsp = find(binnedResponses_global~=0);
y_tsp = binnedResponses_global(binnedResponses_global~=0)';
mov_filtered_spiked = mov_filtered(:,tsp);
mov_filtered_spiked = gpuArray(mov_filtered_spiked);
filters = gpuArray(filters); bias = gpuArray(bias);

% filters = gpuArray(filters);
% bias{filters} = gpyArray(bias);

kx=gpuArray(zeros(nFrontEnds,N));
while(togo==1)

kx = exp(filters*mov_filtered + repmat(bias,[1,N]));


% normalize

kx_sum=sum(kx,1);
kx_sum = kx_sum(tsp);
kx_sum=kx_sum+mu;


% put things on GPU

for ifilter=1:nFrontEnds
[filters(ifilter,:),bias(ifilter)]= fit_filter_cvx(y_tsp,kx(ifilter,:)',tsp,kx_sum,N,mov_filtered_b,dt,mov_filtered_spiked,gather(filters(ifilter,:)),gather(bias(ifilter)));
end


lam = sum(kx,1);
lam=lam+mu;

likelihood = (sum(-lam*(interval/120)) + sum(y_tsp.*log(lam(tsp))))/size(mov_filtered,2);
f_val_prev=f_val;
f_val=-likelihood;
icnt=icnt+1; f_val_log(icnt,1)=f_val;
if((abs(f_val-f_val_prev)/abs(f_val))<tol)
togo=0;
end

end

fitGLM.data_act.kx = kx;
fitGLM.data_act.lam = lam;

for ifilter=1:nFrontEnds
fitGLM.Linear.filter{ifilter}=reshape(gather(filters(ifilter,:)),[filteredStimDim,1]);
fitGLM.Linear.bias{ifilter} = gather(bias(ifilter));
end

f_val=gather(f_val);

fitGLM.mu=mu;
% figure;
% plot(f_val_log);
% title('Objective v/s iterations');
end
