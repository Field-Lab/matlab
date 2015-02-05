function [fitGLM,output] = fitGMLM_afterSTC(binnedResponses,mov,WN_uSq,WNSTA,subunits)


global_vars_GMLM_afterSTC
binnedResponses_global=binnedResponses; 
%%
nSTCs=3;

%% Normalize STA and STC 

WNSTA=WNSTA/norm(WNSTA(:)); % Normalize
for istc=1:nSTCs
    WN_uSq{istc} = WN_uSq{istc}/norm(WN_uSq{istc}(:)); % Normalize
end

%%
mov_filtered = squeeze(convn(mov,WNSTA,'valid'))';
for istc=1:nSTCs

    mov_filtered(istc+1,:)=squeeze(convn(mov,WN_uSq{istc},'valid'))';
end

binnedResponses_global=binnedResponses(end-size(mov_filtered,2)+1:end);

%mov_filtered=[mov_filtered;ones(1,size(mov_filtered,2))]; % Add DC 1?

filteredStimDim=size(mov_filtered,1);
nFrontEnds=4;
temporalCouplingLen=0;

%% Initialize with sub-units to check if the sub-units are optimal weights

% Find weights for sub-units 
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
%initialFilters = (2/nFrontEnds)*(rand(filteredStimDim*nFrontEnds,1)-0.5);
%initialFilters = 2*(rand(filteredStimDim*nFrontEnds,1)-0.5);

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
initialFilters= eye(nFrontEnds,nFrontEnds);
initialFilters=initialFilters(:);

%%
%options = optimoptions(@fminunc,'GradObj','on','Diagnostics','on');

optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','iter',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','off');
 [x,fval,exitflag,output,grad,hessian]  = fminunc(@GMLM_afterSTC,initialFilters,optim_struct);

 change_from_initial = norm(initialFilters-x)

fitGLM.Linear.filter=cell(4,1);
for ifilter=1:nFrontEnds
fitGLM.Linear.filter{ifilter}=x((ifilter-1)*filteredStimDim+1:ifilter*filteredStimDim);
end

FrontEnds=cell(4,1);
for ifilter=1:nFrontEnds
FrontEnds{ifilter}=fitGLM.Linear.filter{ifilter}(1)*WNSTA;
for iSTC=1:nSTCs
    FrontEnds{ifilter} = FrontEnds{ifilter} + fitGLM.Linear.filter{ifilter}(iSTC+1)*WN_uSq{iSTC};
end
end



%% Plot Subunits
figure;
for ifilter=1:nFrontEnds
subplot(3,2,ifilter);
imagesc(subunits{ifilter}(:,:,4));
colormap gray
colorbar
axis image
title('Sub-units');
end

%% Plot STC+STA
figure;
subplot(3,2,1);
imagesc(WNSTA(:,:,4));
colormap gray
colorbar
axis image
title('STA');

for iSTC=1:nSTCs
subplot(3,2,iSTC+1);
imagesc(WN_uSq{iSTC}(:,:,4));
colormap gray
colorbar
axis image
title('STC');
end

%% Plot derived filters
figure;
for ifilter=1:nFrontEnds
subplot(3,2,ifilter);
imagesc(FrontEnds{ifilter}(:,:,4));
colormap gray
colorbar
axis image
title('GMLM filters');
end

end