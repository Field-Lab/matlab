function [fitGLM,output] = fitGMLM_afterSTC_simplified(binnedResponses,mov,filteredStimDim2,nFrontEnds2)


global_vars_GMLM_afterSTC
binnedResponses_global=binnedResponses; 
 filteredStimDim=filteredStimDim2;
 nFrontEnds=nFrontEnds2;

%% Doubt! 
mov_filtered=mov;
% spikeTriggeredSubUnitInput(binnedResponses_global,mov_filtered')

%% Filter with common temporal waveform ??? Better to debug, if we work with spatial part only ? .. valid/same .. doubtful! .. take from SpatialSTC code.. 
%% Fit GMM and find cluster centers
% 
% gm = fitGMM(binnedResponses_global,mov_filtered');
% mu2=gm.mu';
% initialFilters = mu2(:);

%% Initialize with sub-units to check if the sub-units are optimal weights
% 
% %Find weights for sub-units 
% initialFilters = (1/nFrontEnds)*ones(filteredStimDim*nFrontEnds,1);
% nSTCs=3;
% WN_uSq=WN_uSq;
% for iSTC=1:nSTCs
%     WN_uSq{iSTC}=WN_uSq{iSTC}/norm(WN_uSq{iSTC}(:));
% end
% STAnormed=WNSTA/norm(WNSTA(:));
% 
% subUniterr=[];
% icnt=0;
% for isubUnit=1:length(subunits)
%     su=squeeze(subunits{isubUnit});
%       
%     initialFilters((isubUnit-1)*filteredStimDim+1)=sum(su(:).*STAnormed(:));
%     
%     for iSTC=1:nSTCs
%     initialFilters((isubUnit-1)*filteredStimDim+1+iSTC) = sum(su(:).*WN_uSq{iSTC}(:));
%     end
%  
% end
% display('Weights initialized to true subunits');
%% Equal initialization
% initialFilters = (1/nFrontEnds)*ones(filteredStimDim*nFrontEnds,1);

%%  Random initialization
initialFilters = (2/nFrontEnds)*(rand(filteredStimDim*nFrontEnds,1)-0.5);
initialFilters = 2*(rand(filteredStimDim*nFrontEnds,1)-0.5);

%% Random equi-norm initialization
%  initialFilters=zeros(filteredStimDim*nFrontEnds,1);
%  for ifrontend=1:nFrontEnds
%     su=squeeze(subunits{ifrontend});
%       
%     initialFilters((ifrontend-1)*filteredStimDim+1:(ifrontend-1)*filteredStimDim+nFrontEnds)=randn(nFrontEnds,1);
%      initialFilters((ifrontend-1)*filteredStimDim+1:(ifrontend-1)*filteredStimDim+nFrontEnds)= initialFilters((ifrontend-1)*filteredStimDim+1:(ifrontend-1)*filteredStimDim+nFrontEnds)/norm( initialFilters((ifrontend-1)*filteredStimDim+1:(ifrontend-1)*filteredStimDim+nFrontEnds));
%  
%  end
%  display('Random Equinorm initialization')
 %% EYE filters
% initialFilters= eye(nFrontEnds,nFrontEnds);
% initialFilters=initialFilters(:);

%%
%options = optimoptions(@fminunc,'GradObj','on','Diagnostics','on');
%% Add bias in initial filter
initialFilters=[initialFilters;0];
%%
optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','iter',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','off');
 [x,fval,exitflag,output,grad,hessian]  = fminunc(@GMLM_afterSTC,initialFilters,optim_struct);

 change_from_initial = norm(initialFilters-x)/norm(initialFilters)
 mu = x(end)

 fitGLM.Linear.filter=cell(4,1);
for ifilter=1:nFrontEnds
fitGLM.Linear.filter{ifilter}=x((ifilter-1)*filteredStimDim+1:ifilter*filteredStimDim);
end

fitGLM.mu=mu;
end

