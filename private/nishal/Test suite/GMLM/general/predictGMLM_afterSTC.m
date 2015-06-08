function [response,mov_filtered] = predictGMLM_afterSTC(fitGLM,mov,nTrials)

nSTCs=3;

%% Normalize STA and STC 
WNSTA=fitGLM.WNSTA;
WN_uSq = fitGLM.WN_uSq;
WNSTA=WNSTA/norm(WNSTA(:)); % Normalize
for istc=1:nSTCs
    WN_uSq{istc} = WN_uSq{istc}/norm(WN_uSq{istc}(:)); % Normalize
end

%%
mov_filtered = squeeze(convn(mov,WNSTA,'valid'))';
for istc=1:nSTCs

    mov_filtered(istc+1,:)=squeeze(convn(mov,WN_uSq{istc},'valid'))';
end

%mov_filtered=[mov_filtered;ones(1,size(mov_filtered,2))]; % Add DC 1?

filteredStimDim=size(mov_filtered,1);
nFrontEnds=4;

mu=fitGLM.mu;

filters=fitGLM.Linear.filter;
kx=cell(nFrontEnds,1);
grad1=cell(nFrontEnds,1);
for ifilter=1:nFrontEnds
    kx{ifilter} = exp(filters{ifilter}'*mov_filtered);
end

lam=(kx{1}+kx{2}+kx{3}+kx{4}+mu).*0.0083;
noBins = size(lam,2);
response = zeros(nTrials,noBins);
for iTrial=1:nTrials
response(iTrial,:) = (rand(1,noBins)<lam);
end

response=response';
end