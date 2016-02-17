function [new_stas,totalMaskAccept,CellMasks]=clipSTAs_use_mask(mask_loc,vis_ids,stas,cell_params)

nCells=length(stas);
CellNoiseSigmas=[];
CellMaxSTA=[];
CellSNR=[];
CellMasks=cell(nCells,1);
clippedSTAs=cell(nCells,1);
mm = load(mask_loc);

for icell=1:nCells
    icell
    
    close all
dummySTA= stas{icell};
cellID = vis_ids(icell);
mask_rf = mm.mask(:,:,(cell2mat(mm.cells)==cellID));


% figure;
% imagesc(mask_rf)

CellMasks{icell}=mask_rf;

%mask_rf=imdilate(logical(full(sig_stixels)),[1,1,1,1,1;1,1,1,1,1;1,1,1,1,1]);

fastClipSTA=zeros(size(dummySTA));
for itime=1:size(dummySTA,4)
    for idim=1:size(dummySTA,3)
fastClipSTA(:,:,idim,itime)=mask_rf.*dummySTA(:,:,idim,itime);
    end
end
% 
% 
% figure;
% plot(squeeze(sum(sum(sum(fastClipSTA(:,:,:,:),1),2),3)));
% xlabel('STA frame');
% title('Clipped STA summation of pixel values in each frame (Temporal profile)');

cutoff_after=cell_params.STAlen;%input('Cutoff STA after what number? '); %14;
fastClipSTA(:,:,:,cutoff_after:end)=0;

clippedSTAs{icell}=fastClipSTA;

end


newSTAs=cell(nCells,1);
icnt=0;
totalMaskAccept=zeros(size(mask_rf));
totalMaskReject=zeros(size(mask_rf));
for icell=1:nCells
        icnt=icnt+1;
        newSTAs{icnt}=clippedSTAs{icell};
        totalMaskAccept=totalMaskAccept+CellMasks{icell};

end

figure('Color','w');
subplot(1,2,1)
imagesc(totalMaskReject);
title('Rejected RFs');
subplot(1,2,2)
imagesc(totalMaskAccept);
title('Accepted RFs');

new_stas=newSTAs;
%totalMaskAccept=double(totalMaskAccept~=0);
pause(1)
end