function get_su_mask(folder,cellIDs,nSUs,filt_list,thresh,stixSz,maxFrNo)

%% get mask
stas_t = cell(1,1);

icell=1;
for cellID=cellIDs
   data = load(['/Volumes/Lab/Users/bhaishahster/',folder,sprintf('/Cell_%d.mat',cellID)]);
     nSU=nSUs(icell);
    fitGMLM = data.fitGMLM_log{nSU};
  
   
    for ifilt = filt_list
stas = fitGMLM.stas_true{ifilt};
maxfr = (stas(:,:,1,maxFrNo));
maxx =max(maxfr(:));
mask = double(maxfr>thresh*maxx);
if(~exist('totalMask','var'))
    totalMask=mask;
else
    totalMask = totalMask+ mask;
end
    end
    
    
end

totalMask = double(totalMask>0);
iidx = 1:sum(totalMask(:));
totalMask(totalMask>0) = iidx;

totalMask = repelem(totalMask,stixSz,stixSz);
%% write mask
dlmwrite([['/Volumes/Lab/Users/bhaishahster/',folder,sprintf('/mask_sus_%d,%d.mat',filt_list(1),filt_list(2))] '.txt'], totalMask, 'delimiter', '\t', 'newline', 'pc'); % if this errors which save path
savedMap = dlmread([['/Volumes/Lab/Users/bhaishahster/',folder,sprintf('/mask_sus_%d,%d.mat',filt_list(1),filt_list(2))] '.txt']);

%% ------------------------------- Display mask ----------------------------------------------
figure
imagesc(savedMap)
title('Stixels to be modulated')
axis equal
end