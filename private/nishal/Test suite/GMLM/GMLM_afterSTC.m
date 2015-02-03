function [fval,Grad] = GMLM_afterSTC(x)

global_vars_GMLM_afterSTC


ers=cell(nFrontEnds,1);
for ifilter=1:nFrontEnds
filters{ifilter}=x((ifilter-1)*filteredStimDim+1:ifilter*filteredStimDim);
end

kx=cell(nFrontEnds,1);
grad1=cell(nFrontEnds,1);
for ifilter=1:nFrontEnds
    kx{ifilter} = exp(filters{ifilter}'*mov_filtered);
grad1{ifilter}=-sum(repmat(kx{ifilter},[filteredStimDim,1]).*mov_filtered,2)*0.0083; %(dt = 0.0083)
end


tsp = find(binnedResponses_global~=0);
kx_sum=0*kx{1}(tsp);
for ifilter=1:nFrontEnds
kx_sum = kx_sum + kx{ifilter}(tsp);
end

grad2=cell(nFrontEnds,1);
mov_filtered_spiked = mov_filtered(:,tsp);
for ifilter=1:nFrontEnds
grad2{ifilter}= sum(repmat((kx{ifilter}(tsp)./kx_sum),[filteredStimDim,1]).*mov_filtered_spiked,2);
end

grad_K =zeros(filteredStimDim*nFrontEnds,1);
for ifilter=1:nFrontEnds
grad_K((ifilter-1)*filteredStimDim+1:ifilter*filteredStimDim)=grad1{ifilter}+grad2{ifilter};
end
Grad=-grad_K;

grad_norm= norm(Grad);


lam=(kx{1}+kx{2}+kx{3}+kx{3});
likelihood = sum(-lam*0.0083) + sum(log(lam(tsp)));
fval=-likelihood;

end