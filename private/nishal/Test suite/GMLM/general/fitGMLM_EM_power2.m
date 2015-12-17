function [fitGLM,f_val_log,f_val_test_log,su_compare_log ] = fitGMLM_EM_power2(binnedResponses,mov,filteredStimDim2,nFrontEnds2,interval,gamma,varargin)

%% Works only for gaussian stimuli and no bias inside exponential 

%% Data
global_vars_GMLM_afterSTC
 binnedResponses_global=binnedResponses; 
 filteredStimDim=filteredStimDim2;
 nFrontEnds=nFrontEnds2;
 
mov_filtered=mov;

%% parse input
 p =inputParser;
% specify list of optional parameters
%initialFilters = (2/nFrontEnds)*(rand(filteredStimDim*nFrontEnds,1)-0.5);
initialFilt = 2*(rand(filteredStimDim*nFrontEnds,1)-0.5);
p.addParamValue('initialFilters',initialFilt);

testStruct_def.test=0;
p.addParamValue('testStruct',testStruct_def);

comp_true_SU_str.compare=0;
p.addParamValue('compare_true_SU',comp_true_SU_str);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% assign variables from parsing.
initialFilters = params.initialFilters;
testStruct = params.testStruct;
compare_true_SU=params.compare_true_SU;


%% Initialization.
%  Random initialization
% initialFilters = (2/nFrontEnds)*(rand(filteredStimDim*nFrontEnds,1)-0.5);
% initialFilters = 2*(rand(filteredStimDim*nFrontEnds,1)-0.5);
 
% Add bias in initial filter
%%%%%% initialFilters=[initialFilters;0];
 

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
mov_len=size(mov_filtered,2);

kx = gpuArray(zeros(nFrontEnds,mov_len));
kkx = gpuArray(zeros(nFrontEnds,mov_len));
binnedResponses_global = gpuArray(binnedResponses_global);

f_val_log=gpuArray([]); icnt=0;f_val_test_log=[];su_compare_log=[];
iter=0;
while(togo==1 & iter<=max_iter)
iter=iter+1;
%    kx=cell(nFrontEnds,1);
%kkx = cell(nFrontEnds,1);

if(testStruct.test~=0)
f_val_test_log = [f_val_test_log;test_data(testStruct,filters,gamma,interval,mu)];
end

if(compare_true_SU.compare~=0)
su_compare_log = [su_compare_log;compare_SU(filters,compare_true_SU)];
end

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

lam=sum(kx,1);
lam=lam+mu;

likelihood = (sum(-lam*(interval/120)) + sum(y_tsp.*log(lam(tsp))))/size(mov_filtered,2);
f_val_prev=f_val;
f_val=-likelihood;

icnt=icnt+1; f_val_log(icnt,1)=f_val;
if((abs(f_val-f_val_prev)/abs(f_val))<tol)
togo=0;
end

end
if(iter==max_iter)
display('Stopping because max iter reached');    
end

phase2=0;
if(phase2==1)
    
k_est=filters;
  
    su_mask=double((-k_est)==repmat(max((-k_est),[],1),[nFrontEnds,1])); display('OFF cell assumed');
    
    
    togo=1; iter=0;
    kx=gpuArray(zeros(nFrontEnds,size(mov_filtered,2)));
    kkx = kx;
    while(togo==1 & iter<=max_iter)
        iter=iter+1;

filters = filters.*su_mask;
kkx= filters*mov_filtered;
kx = (kkx.*(kkx>0)).^gamma;

if(testStruct.test~=0)
f_val_test_log = [f_val_test_log;test_data(testStruct,filters,gamma,interval,mu)];
end

if(compare_true_SU.compare~=0)
su_compare_log = [su_compare_log;compare_SU(filters,compare_true_SU)];
end
% normalize
tsp = find(binnedResponses_global~=0);
y_tsp = binnedResponses_global(binnedResponses_global~=0)';
kx_total_sum = sum(kx,1);

kx_sum = kx_total_sum(tsp);
kx_sum=kx_sum+mu;

mov_filtered_spiked = mov_filtered(:,tsp);

for ifilter=1:nFrontEnds
    kkxtsp = kkx(ifilter,tsp);
    mask = logical(su_mask(ifilter,:)');
    if(sum(mask)==0)
    continue;
    end
    mov_fs = mov_filtered_spiked(mask,:);
    kkxtsp_gam=y_tsp.*(kkxtsp>0).*((kkxtsp).^(gamma-1))./kx_sum;
sta_f = gamma*mov_fs*kkxtsp_gam'/size(mov_fs,2);
% filters{ifilter} = (Sigma\sta_f{ifilter})/N;
optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','off',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','on');

[x,fval,exitflag,output,grad,hessian]  = fminunc(@(x)GMLM_quad_decomposed_fcn(x,mov_filtered(mask,:),sta_f,size(mov_filtered,2),interval/120,gamma),gather(filters(ifilter,mask)'),optim_struct);
filters(ifilter,mask) = gpuArray(x(1:end));
end

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
display('Phase 2 Stopping because max iter reached');    
end

    
end
filters = gather(filters);
for ifilter= 1:nFrontEnds
fitGLM.Linear.filter{ifilter}=filters(ifilter,:)';
end

fitGLM.mu=mu;
fitGLM.data_act.kx=gather(kx);
fitGLM.data_act.lam=gather(lam);
% 
%  figure;
%  plot(f_val_log);
%  title('Objective v/s iterations');
end

function f_val = test_data(testStruct,filters,gamma,interval,mu)

binnedResponses  =testStruct.binnedResponses;
mov_filtered = testStruct.mov;



    kkx= filters*mov_filtered;
    kx = (kkx.*(kkx>0)).^gamma;
    tsp = find(binnedResponses~=0);

lam=sum(kx,1);
lam=lam+mu;

y_tsp =binnedResponses(binnedResponses~=0)';
likelihood = (sum(-lam*(interval/120)) + sum(y_tsp.*log(lam(tsp))))/size(mov_filtered,2);

f_val=-likelihood;
end

function dist = compare_SU(filters,comp_true_SU_str)
filters= filters';
su = comp_true_SU_str.su; su=su';
proj = comp_true_SU_str.proj;

filters = proj*filters;
su = proj*su;
max_correlation = max(abs(corr(gather(filters),gather(su))),[],2);
dist = sum(max_correlation);

end

