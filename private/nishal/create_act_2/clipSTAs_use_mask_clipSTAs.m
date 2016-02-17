function [new_stas,totalMaskAccept,CellMasks]=clipSTAs_use_mask_clipSTAs(mask_loc,vis_ids);

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
mask_rf = mm.mask(:,:,(mm.cells==vis_ids))


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

CellNoiseSigmas(icell)=params.noise_sigmas(1); % TODO nishal
CellMaxSTA(icell)=(max(abs(fastClipSTA(:))));

clippedSTAs{icell}=fastClipSTA;

end

CellSNR=CellMaxSTA./CellNoiseSigmas;
figure;
hist(CellSNR,20);

threshold=-1%input('Select SNR threshold to eliminate cells?');

if(threshold>=0)
cellChoose=(CellSNR>threshold);
else
cellChoose=ones(size(CellSNR));    
end

newSTAs=cell(sum(cellChoose),1);
icnt=0;
totalMaskAccept=zeros(size(mask_rf));
totalMaskReject=zeros(size(mask_rf));
for icell=1:nCells
    if(cellChoose(icell))
        icnt=icnt+1;
        newSTAs{icnt}=clippedSTAs{icell};
        totalMaskAccept=totalMaskAccept+CellMasks{icell};
    else
        display(sprintf('Cell %d discarded',icell));
        totalMaskReject=totalMaskReject+CellMasks{icell};
    end
    
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