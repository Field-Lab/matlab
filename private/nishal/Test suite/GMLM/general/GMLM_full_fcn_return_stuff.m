function [fval,Grad,lam] = GMLM_full_fcn_return_stuff(x)

global_vars_GMLM_afterSTC
dt=1/1200;
%% spontaneous activity
mu=x(end);
mu=0; %% mu stuff!
%display('ALERT : NoMu');

%% Stimulus filters
x=x(1:end-1);

filters=cell(nFrontEnds,1);
for ifilter=1:nFrontEnds
filters{ifilter}=x((ifilter-1)*filteredStimDim+1:ifilter*filteredStimDim);
end
histBasLen=size(hInp,1);
histBas = x(nFrontEnds*filteredStimDim+1:nFrontEnds*filteredStimDim + histBasLen);

kh = exp(histBas'*hInp);

% grad 1 terms
kx=cell(nFrontEnds,1);
grad1=cell(nFrontEnds,1);
for ifilter=1:nFrontEnds
    kx{ifilter} = exp(filters{ifilter}'*mov_filtered);
   
grad1{ifilter}=-(mov_filtered*(kx{ifilter}.*kh)')*dt; %(dt = 0.0083)
end
grad1_mu=-sum(kh)*dt; %(dt = 0.0083)

% grad 2 terms
tsp = find(binnedResponses_global~=0);
y_tsp = binnedResponses_global(binnedResponses_global~=0)';

kx_sum=0*kx{1}(tsp);
for ifilter=1:nFrontEnds
kx_sum = kx_sum + kx{ifilter}(tsp);
end
kx_sum=kx_sum+mu;

grad2=cell(nFrontEnds,1);
mov_filtered_spiked = mov_filtered(:,tsp);
for ifilter=1:nFrontEnds 
grad2{ifilter}= mov_filtered_spiked*( y_tsp.*kx{ifilter}(tsp)./kx_sum)';
end
grad2_mu = sum((y_tsp./kx_sum),2);


lam=kx{1};
for ik=2:length(kx)
lam=lam+kx{ik};
end
lam=lam+mu;
lam = lam.*kh;


grad1_h=-hInp*lam'*dt; %(dt = 0.0083)
grad2_h = hInp*binnedResponses_global;


grad_K =zeros(filteredStimDim*nFrontEnds+1+histBasLen,1);
for ifilter=1:nFrontEnds
grad_K((ifilter-1)*filteredStimDim+1:ifilter*filteredStimDim)=grad1{ifilter}+grad2{ifilter};
end
% grad_K=grad_K*0; %display('Filters Not changed');
grad_K(end) = 0; %grad1_mu+grad2_mu;
grad_K(nFrontEnds*filteredStimDim+1:nFrontEnds*filteredStimDim + histBasLen)=grad1_h + grad2_h;

Grad=-grad_K/size(mov_filtered,2);

grad_norm= norm(Grad);

likelihood = sum(-lam*dt) + sum(y_tsp.*log(lam(tsp)));
fval=-likelihood/size(mov_filtered,2);

end