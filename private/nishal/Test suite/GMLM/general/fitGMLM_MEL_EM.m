function [fitGLM,f_val] = fitGMLM_MEL_EM(binnedResponses,mov,filteredStimDim2,nFrontEnds2,interval)

%% Works only for gaussian stimuli and no bias inside exponential 
 
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
for ifilter=1:nFrontEnds
filters{ifilter}=x((ifilter-1)*filteredStimDim+1:ifilter*filteredStimDim);
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
Sigma = mov_filtered*mov_filtered'/size(mov_filtered,2)
Sigmainv=inv(Sigma);

f_val_log=[]; icnt=0;
while(togo==1)
kx=cell(nFrontEnds,1);

for ifilter=1:nFrontEnds
    kx{ifilter} = exp(filters{ifilter}'*mov_filtered);
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
sta_f{ifilter}= sum(repmat((y_tsp.*kx{ifilter}(tsp)./kx_sum),[filteredStimDim,1]).*mov_filtered_spiked,2);
C = norm(Sigma\sta_f{ifilter})/N;
filters{ifilter}= find_alpha(C,Sigma) * (Sigma\sta_f{ifilter})/norm(Sigma\sta_f{ifilter});
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


if(phase2==1)
togo=1;
display('Phase 2');
    while(togo==1)
    
    mu  = sum(binnedResponses)/( size(mov_filtered,2)*interval/120);
    for ifilter=1:nFrontEnds
    mu = mu - exp((0.5)*filters{ifilter}'*Sigma*filters{ifilter});
    end
    mu=max(mu,0);
    
kx=cell(nFrontEnds,1);

for ifilter=1:nFrontEnds
    kx{ifilter} = exp(filters{ifilter}'*mov_filtered);
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
sta_f{ifilter}= sum(repmat((y_tsp.*kx{ifilter}(tsp)./kx_sum),[filteredStimDim,1]).*mov_filtered_spiked,2);
C = norm(Sigma\sta_f{ifilter})/N;
filters{ifilter}= find_alpha(C,Sigma) * (Sigma\sta_f{ifilter})/norm(Sigma\sta_f{ifilter});
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

phase3=0;
if(phase3==1)
    initialFilters=zeros(nFrontEnds*filteredStimDim + 1,1);
    initialFilters(end)=mu;
    for ifilter=1:nFrontEnds
    initialFilters((ifilter-1)*filteredStimDim +1 : ifilter*filteredStimDim,1) = filters{ifilter};
    end
    
    
optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','iter',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','off');
 [x,fval,exitflag,output,grad,hessian]  = fminunc(@GMLM_afterSTC,initialFilters,optim_struct);
 mu = x(end)

filters=cell(nFrontEnds,1);
for ifilter=1:nFrontEnds
filters{ifilter}=x((ifilter-1)*filteredStimDim+1:ifilter*filteredStimDim);
end

f_val_log = [f_val_log,fval];
end


fitGLM.Linear.filter=filters;
fitGLM.mu=mu;
figure;
plot(f_val_log);
title('Objective v/s iterations');
end

function a= find_alpha(C,Sigma)
al=-1; au=C;
tol=1e-5;
togo=1;
while(togo==1)
a=(al+au)/2;

v=a*exp((0.5)*(a^2*Sigma(1,1))) - C ;

if(v>0)
    au=a;
else
    al=a;
end
if(abs(v)<tol)
togo=0;
end

end
end