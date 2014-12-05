function [new_stas,totalMaskAccept]=clipSTAs(stas,cell_params)

nCells=length(stas);
CellNoiseSigmas=[];
CellMaxSTA=[];
CellSNR=[];
CellMasks=cell(nCells,1);
clippedSTAs=cell(nCells,1);

for icell=1:nCells
    icell
    
    close all
dummySTA= stas{icell};


[sig_stixels, params, rf_strength_out] = significant_stixels(dummySTA,'select','thresh','thresh',2);
% figure;
% imagesc(sig_stixels)    



% Find continous region
max_val=double(max(rf_strength_out(sig_stixels)));
[row,col]=find(double(rf_strength_out).*double(sig_stixels)==max_val)
list=[row,col];

[row,col]=find(sig_stixels==1);
sig_list=[row,col];


tostop=0;
while tostop==0
tostop=1;
for istix=1:size(sig_list,1)
    if ~(sum((sig_list(istix,1)==list(:,1)) .* (sig_list(istix,2)==list(:,2)))==1)
        for ilist_elem=1:size(list,1)
            if(abs(sig_list(istix,1)-list(ilist_elem,1)) + abs(sig_list(istix,2)-list(ilist_elem,2))<=1)
            list=[list;sig_list(istix,1),sig_list(istix,2)];
            tostop=0;
%             display('Element added')
            a=sig_list(istix,:);
            b=list(ilist_elem,:);
            break
            end
        end
    end
end
end

mask_rf=zeros(size(sig_stixels));
for ielem=1:size(list,1)
mask_rf(list(ielem,1),list(ielem,2))=1;
end

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

threshold=input('Select SNR threshold to eliminate cells?');


cellChoose=(CellSNR>threshold);

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

figure;
subplot(1,2,1)
imagesc(totalMaskReject);
title('Rejected RFs');
subplot(1,2,2)
imagesc(totalMaskAccept);
title('Accepted RFs');

new_stas=newSTAs;
totalMaskAccept=double(totalMaskAccept~=0);
pause(1)
end