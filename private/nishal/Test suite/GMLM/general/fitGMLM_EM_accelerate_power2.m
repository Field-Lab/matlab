function [fitGLM,f_val] = fitGMLM_EM_accelerate_power2(binnedResponses,mov,filteredStimDim2,nFrontEnds2,interval,gamma)

%% Acceleration rate
alpha = 0.2;

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

filters=zeros(nFrontEnds,filteredStimDim);
for ifilter=1:nFrontEnds
filters(ifilter,:)=x((ifilter-1)*filteredStimDim+1:ifilter*filteredStimDim);
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
max_iter=1000;
% make 0 mean. 

mov_filtered= gpuArray(mov_filtered);
filters = gpuArray(filters);
filters_prev = filters;
mov_len=size(mov_filtered,2);

kx = gpuArray(zeros(nFrontEnds,mov_len));
kkx = gpuArray(zeros(nFrontEnds,mov_len));
binnedResponses_global = gpuArray(binnedResponses_global);

f_val_log=gpuArray([]); icnt=0;
iter=0;
while(togo==1 & iter<=max_iter)
iter=iter+1
%    kx=cell(nFrontEnds,1);
%kkx = cell(nFrontEnds,1);

    kkx= filters*mov_filtered;
    kx = (kkx.*(kkx>0)).^gamma;


% normalize
tsp = find(binnedResponses_global~=0);
y_tsp = binnedResponses_global(binnedResponses_global~=0)';
kx_total_sum = sum(kx,1);

kx_sum = kx_total_sum(tsp);
kx_sum=kx_sum+mu;



mov_filtered_spiked = mov_filtered(:,tsp);


for ifilter=1:nFrontEnds
    kkxtsp = kkx(ifilter,tsp);
    kkxtsp_gam=y_tsp.*(kkxtsp>0).*((kkxtsp).^(gamma-1))./kx_sum;
sta_f = gamma*mov_filtered_spiked*kkxtsp_gam'/size(mov_filtered,2);
% filters{ifilter} = (Sigma\sta_f{ifilter})/N;
optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','off',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','on');

[x,fval,exitflag,output,grad,hessian]  = fminunc(@(x)GMLM_quad_decomposed_fcn(x,mov_filtered,sta_f,size(mov_filtered,2),interval/120,gamma),gather(filters(ifilter,:)'),optim_struct);
filters(ifilter,:) = gpuArray(x(1:end));
end

filters = (1-alpha)*filters + alpha*filters_prev;
filters_prev = filters;


lam=sum(kx,1);
lam=lam+mu;

likelihood = (sum(-lam*(interval/120)) + sum(y_tsp.*log(lam(tsp))))/size(mov_filtered,2);
f_val_prev=f_val;
f_val=-likelihood
icnt=icnt+1; f_val_log(icnt,1)=f_val;
if((abs(f_val-f_val_prev)/abs(f_val))<tol)
togo=0;
end

end
if(iter==max_iter)
display('Stopping because max iter reached');    
end


    
filters = gather(filters);
for ifilter= 1:nFrontEnds
fitGLM.Linear.filter{ifilter}=filters(ifilter,:)';
end

fitGLM.mu=mu;
fitGLM.data_act.kx=gather(kx);
fitGLM.data_act.lam=gather(lam);
% 
 figure;
 plot(f_val_log);
 title('Objective v/s iterations');
end
