function [fval,Grad] = GMLM_afterSTC(x)

global_vars_GMLM_afterSTC

%% spontaneous activity
mu=x(end);
mu=0; %% mu stuff!

%% Stimulus filters
x=x(1:end-1);

filters=cell(nFrontEnds,1);
for ifilter=1:nFrontEnds
filters{ifilter}=x((ifilter-1)*filteredStimDim+1:ifilter*filteredStimDim);
end

% grad 1 terms
kx=cell(nFrontEnds,1);
grad1=cell(nFrontEnds,1);
for ifilter=1:nFrontEnds
    kx{ifilter} = exp(filters{ifilter}'*mov_filtered);
grad1{ifilter}=-sum(repmat(kx{ifilter},[filteredStimDim,1]).*mov_filtered,2)*0.0083; %(dt = 0.0083)
end
grad1_mu=-(size(mov_filtered,2))*0.0083; %(dt = 0.0083)

% grad 2 terms
tsp = find(binnedResponses_global~=0);
kx_sum=0*kx{1}(tsp);
for ifilter=1:nFrontEnds
kx_sum = kx_sum + kx{ifilter}(tsp);
end
kx_sum=kx_sum+mu;

grad2=cell(nFrontEnds,1);
mov_filtered_spiked = mov_filtered(:,tsp);
for ifilter=1:nFrontEnds
grad2{ifilter}= sum(repmat((kx{ifilter}(tsp)./kx_sum),[filteredStimDim,1]).*mov_filtered_spiked,2);
end
grad2_mu = sum((1./kx_sum),2);



grad_K =zeros(filteredStimDim*nFrontEnds+1,1);
for ifilter=1:nFrontEnds
grad_K((ifilter-1)*filteredStimDim+1:ifilter*filteredStimDim)=grad1{ifilter}+grad2{ifilter};
end
grad_K(nFrontEnds*filteredStimDim+1) = grad1_mu+grad2_mu;

Grad=-grad_K;

grad_norm= norm(Grad);


lam=(kx{1}+kx{2}+kx{3}+kx{4}+mu);
likelihood = sum(-lam*0.0083) + sum(log(lam(tsp)));
fval=-likelihood;

end