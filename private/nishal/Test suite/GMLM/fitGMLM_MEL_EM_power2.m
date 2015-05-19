function [fitGLM,f_val] = fitGMLM_MEL_EM_power2(binnedResponses,mov,filteredStimDim2,nFrontEnds2,interval,gamma)

%% Works only for gaussian stimuli and no bias inside exponential 
 
%% Data
global_vars_GMLM_afterSTC
 binnedResponses_global=binnedResponses; 
 filteredStimDim=filteredStimDim2;
 nFrontEnds=nFrontEnds2;
 
mov_filtered=mov;


%% Initialization.
%  Random initialization
initialFilters = (2/nFrontEnds)*(rand(filteredStimDim*nFrontEnds,1)-0.5);
initialFilters = 2*(rand(filteredStimDim*nFrontEnds,1)-0.5);
 % Add bias in initial filter
% initialFilters=[initialFilters;0];
 

% Extract filters 
x=initialFilters;

filters=cell(nFrontEnds,1);
for ifilter=1:nFrontEnds
filters{ifilter}=x((ifilter-1)*filteredStimDim+1:ifilter*filteredStimDim);
end
%mu=sum(binnedResponses)/(((nFrontEnds +1)^2)* size(mov_filtered,2)*interval/120);
mu=0.1;
%% E-M MEL . 
N =size(mov_filtered,2)*interval/120;
% calculate alphas 
% SU outputs
togo=1;
f_val=Inf;
tol=1e-5;
% make 0 mean. 
Sigma = mov_filtered*mov_filtered'/size(mov_filtered,2)
Sigmainv=inv(Sigma);

f_val_log=[]; icnt=0;
while(togo==1)
kx=cell(nFrontEnds,1);
kkx = cell(nFrontEnds,1);

for ifilter=1:nFrontEnds
    kkx{ifilter}= filters{ifilter}'*mov_filtered;
    kx{ifilter} = (kkx{ifilter}.*(kkx{ifilter}>0)).^gamma;
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
    kkxtsp = kkx{ifilter}(tsp);
    kkxtsp_gam=y_tsp.*(kkxtsp>0).*((kkxtsp).^(gamma-1))./kx_sum;
sta_f{ifilter} = gamma*mov_filtered_spiked*kkxtsp_gam';
filters{ifilter} = (Sigma\sta_f{ifilter})/N;
end


lam=kx{1};
for ik=2:length(kx)
lam=lam+kx{ik};
end
lam=lam+mu;

likelihood = (sum(-lam*(interval/120)) + sum(y_tsp.*log(lam(tsp))))/size(mov_filtered,2);
f_val_prev=f_val;
f_val=-likelihood
icnt=icnt+1; f_val_log(icnt)=f_val;
if((abs(f_val-f_val_prev)/abs(f_val))<tol)
togo=0;
end

end

phase2=1;
if(phase2==1)

    nc=filteredStimDim; nFilters = nFrontEnds;
  k_est = zeros(nFilters,nc);
    for isu=1:nFilters
    k_est(isu,:)=filters{isu}(1:nc);
    end
    
  
    su_mask=double((-k_est)==repmat(max((-k_est)),[nFilters,1])); display('OFF cell assumed');
    
    
    togo=1;
    while(togo==1)
kx=cell(nFrontEnds,1);
kkx = cell(nFrontEnds,1);

for ifilter=1:nFrontEnds
   filters{ifilter} = filters{ifilter}.*su_mask(ifilter,:)';
    kkx{ifilter}= filters{ifilter}'*mov_filtered;
    kx{ifilter} = (kkx{ifilter}.*(kkx{ifilter}>0)).^gamma;
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
    kkxtsp = kkx{ifilter}(tsp);
    mask=logical(su_mask(ifilter,:)');
    mov_fs = mov_filtered_spiked(mask,:);
    
    kkxtsp_gam=y_tsp.*(kkxtsp>0).*((kkxtsp).^(gamma-1))./kx_sum;
sta_f{ifilter} = gamma*mov_fs*kkxtsp_gam';
filters{ifilter}(mask) = (Sigma(mask,mask)\sta_f{ifilter})/N;
end


lam=kx{1};
for ik=2:length(kx)
lam=lam+kx{ik};
end
lam=lam+mu;

likelihood = (sum(-lam*(interval/120)) + sum(y_tsp.*log(lam(tsp))))/size(mov_filtered,2);
f_val_prev=f_val;
f_val=-likelihood
icnt=icnt+1; f_val_log(icnt)=f_val;
if((abs(f_val-f_val_prev)/abs(f_val))<tol)
togo=0;
end

    end
    
    
end
fitGLM.Linear.filter=filters;
fitGLM.mu=mu;
% figure;
% plot(f_val_log);
% title('Objective v/s iterations');
end
