%% Sub-unit and STC reconstruction error. Spatial Only
% Assumes STC is rank 1
% 
% uSq=uSq150min;
% for iSTC=1:4
%     uSq{iSTC}=uSq{iSTC}(:,:,4);
%     uSq{iSTC}=uSq{iSTC}/norm(uSq{iSTC}(:));
% end
% 
% figure;
% subUniterr=[];
% icnt=0;
% for isubUnit=1:nSubunits
%     su=squeeze(subunits{isubUnit}(:,:,1,4));
%     recons=0*su;
%     for iSTC=1:4
%     recons=recons+sum(su(:).*uSq{iSTC}(:))*uSq{iSTC};
%     end
%     subUniterr(isubUnit)=norm(su(:)-recons(:))/norm(su(:));
%     
%     icnt=icnt+1;
%     subplot(nSubunits,2,icnt);
%     imagesc(su);
%     colormap gray;colorbar;
%     title(sprintf('Actual Subunit %d',isubUnit));
%     axis image
%     
%     icnt=icnt+1;
%     subplot(nSubunits,2,icnt);
%     imagesc(recons);
%     colormap gray;colorbar;
%     title(sprintf('Reconstruction Subunit:%d',isubUnit));
%     axis image
% end
    
%% Sub-unit and STC reconstruction error. Spatio-Temporal reconstruction
% Does not assume rank 1 STC.
nSTCs=3;
uSq=uSq;
for iSTC=1:nSTCs
    uSq{iSTC}=uSq{iSTC}/norm(uSq{iSTC}(:));
end

subUniterr=[];
icnt=0;
for isubUnit=1:nSubunits
    su=squeeze(subunits{isubUnit});
    recons=0*su;
    for iSTC=1:nSTCs
    recons=recons+sum(su(:).*uSq{iSTC}(:))*uSq{iSTC};
    end
    subUniterr(isubUnit)=norm(su(:)-recons(:))/norm(su(:));
    
    figure;
    for itime=1:30
        itime
    icnt=icnt+1;
    subplot(1,2,1);
    imagesc(su(:,:,4));
    colormap gray;colorbar;
    title(sprintf('Actual Subunit %d',isubUnit));
    axis image
    
    icnt=icnt+1;
    subplot(1,2,2);
    imagesc(recons(:,:,4));
    colormap gray;colorbar;
    title(sprintf('Reconstruction Subunit:%d',isubUnit));
    axis image
    pause(1/120);
    end

end

% Problem with mean and other stuff???
