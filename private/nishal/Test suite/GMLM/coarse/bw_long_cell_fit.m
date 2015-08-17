% Load data and fitGMLM on data.. 

cell_list = [467,497,1292,1321,1471,1771,1786,2211,2506,2851,3361,3676,3843,4336,4367,4576,5086,5161,5386,5851,5986,6196,6256,7066,7292,7606]

for icell=2:length(cell_list)
    close all
    
load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2012_08_09_3/longBW/cell_%d.mat',cell_list(icell)));
interval=1;
filteredStimDim = size(maskedMovdd_sliced,1);
fitGMLM_log=cell(15,1);
f_val = zeros(15,1);
for nSU =1:15
[fitGMLM,f_val] = fitGMLM_MEL_EM_bias(binnedSpikeResponses',maskedMovdd_sliced*1000,filteredStimDim,nSU,interval); 
fitGMLM_log{nSU} = fitGMLM;
f_val(nSU) = f_val;

pause(0.2)
end
save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2012_08_09_3/longBW/fitcell_%d.mat',cell_list(icell)),'fitGMLM_log','f_val');
end

%% Look at filters ..

icell=2;
load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2012_08_09_3/longBW/fitcell_%d.mat',cell_list(icell)));
load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2012_08_09_3/longBW/cell_%d.mat',cell_list(icell)));
nSU = 5;

szz =[1,1;
    1,2;
    2,2;
    2,2;
    2,3;
    2,3;
    3,3;
    3,3;
    3,3;
    3,4;
    3,4;
    3,4;
    4,4;
    4,4;
    4,4];

 fitGMLM = fitGMLM_log{nSU};
 mask2=totalMaskSliced;
sta_dim1 = size(mask2,1);
sta_dim2 = size(mask2,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask2));



idx = [1:size(STA,1)]';
frame = repmat(idx,[1,size(STA,2)]);
dim1pix = frame (logical(totalMaskSliced));
minpix1= min(dim1pix);
maxpix1 = max(dim1pix);

idx = [1:size(STA,2)];
frame = repmat(idx,[size(STA,1),1]);
dim2pix = frame (logical(totalMaskSliced));
minpix2= min(dim2pix);
maxpix2 = max(dim2pix);


h=figure;
for ifilt=1:nSU
subplot(szz(nSU,1),szz(nSU,2),ifilt)
u_spatial = reshape_vector(fitGMLM.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);
imagesc(u_spatial(minpix1:maxpix1,minpix2:maxpix2));
colormap gray
axis image
set(gca,'xTick',[]);
set(gca,'yTick',[]);
title(sprintf('SU # %d',ifilt));
end

%print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/Exp_ASM_MEL_EM_filters SU%d.pdf',nSU));
