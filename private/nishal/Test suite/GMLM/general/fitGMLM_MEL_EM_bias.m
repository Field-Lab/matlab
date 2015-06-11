function [fitGLM,f_val] = fitGMLM_MEL_EM_bias(binnedResponses,mov,filteredStimDim2,nFrontEnds2,interval)

%% Works only for gaussian stimuli and there is bias inside exponential 
 
%% Data
global_vars_GMLM_afterSTC
 binnedResponses_global=binnedResponses; 
 filteredStimDim=filteredStimDim2;
 nFrontEnds=nFrontEnds2;
 
mov_filtered=mov;
phase2=0;

%% Initialization.
%  Random initialization
initialFilters = (2/nFrontEnds)*(rand(filteredStimDim*nFrontEnds,1)-0.5);
initialFilters = 2*(rand(filteredStimDim*nFrontEnds,1)-0.5);
 % Add bias in initial filter
% initialFilters=[initialFilters;0];
 

% Extract filters 
x=initialFilters;

filters=cell(nFrontEnds,1);
bias=cell(nFrontEnds,1);
for ifilter=1:nFrontEnds
filters{ifilter}=x((ifilter-1)*filteredStimDim+1:ifilter*filteredStimDim);
bias{ifilter}=2*(rand(1,1)-0.5);
end
%mu=sum(binnedResponses)/(((nFrontEnds +1)^2)* size(mov_filtered,2)*interval/120);
mu=0;
%% E-M MEL . 
N =size(mov_filtered,2)*interval/120;
% calculate alphas 
% SU outputs
togo=1;
f_val=Inf;
tol=1e-5;
% make 0 mean. 
Sigma = mov_filtered*mov_filtered'/size(mov_filtered,2);
Sigmainv=inv(Sigma);

f_val_log=[]; icnt=0;
while(togo==1)
kx=cell(nFrontEnds,1);

for ifilter=1:nFrontEnds
    kx{ifilter} = exp(filters{ifilter}'*mov_filtered + bias{ifilter});
end

% normalize
tsp = find(binnedResponses_global~=0);
y_tsp = binnedResponses_global(binnedResponses_global~=0)';

kx_sum=0*kx{1}(tsp);
for ifilter=1:nFrontEnds
kx_sum = kx_sum + kx{ifilter}(tsp);
end
kx_sum=kx_sum+mu;


mov_filtered_spiked = mov_filtered(:,tsp);
for ifilter=1:nFrontEnds
xx = y_tsp.*kx{ifilter}(tsp)./kx_sum;
sta_f{ifilter}= mov_filtered_spiked*xx';
mean_ass_f{ifilter}= sum(xx,2);
filters{ifilter}=  (Sigma\sta_f{ifilter})/mean_ass_f{ifilter};
bias{ifilter} = log((1/N)*mean_ass_f{ifilter}) - (0.5)*filters{ifilter}'*Sigma*filters{ifilter};
end


lam=kx{1};
for ik=2:length(kx)
lam=lam+kx{ik};
end
lam=lam+mu;

likelihood = (sum(-lam*(interval/120)) + sum(y_tsp.*log(lam(tsp))))/size(mov_filtered,2);
f_val_prev=f_val;
f_val=-likelihood;
icnt=icnt+1; f_val_log(icnt)=f_val;
if((abs(f_val-f_val_prev)/abs(f_val))<tol)
togo=0;
end

end

fitGLM.data_act.kx = kx;
fitGLM.data_act.lam = lam;

fitGLM.Linear.filter=filters;
fitGLM.Linear.bias = bias;
fitGLM.mu=mu;
figure;
plot(f_val_log);
title('Objective v/s iterations');
end
