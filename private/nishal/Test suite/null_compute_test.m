%% Calculate null stimulus -in BW
Filtdim1=size(STA,1);
Filtdim2=size(STA,2);

% Select movie for null stim 
movieLen=120*10;
mov_gen2=zeros(Filtdim1,Filtdim2,movieLen);


generate_null_script

 
%% Generate response to null stimulus
stas=cell(1,1);
stas{1}=zeros(Filtdim1,Filtdim2,1,Filtlen);

for itime=1:Filtlen
    if(sta_test==1)
    stas{1}(:,:,1,itime) = squeeze(STA(:,:,itime));
    end
    
    if(sta_test==2)
    stas{1}(:,:,1,itime) = squeeze(sum(full_fit(:,:,:,end-itime+1),3))';
    end
     if(sta_test==3)
    stas{1}(:,:,1,itime) = squeeze(clippedSTA(:,:,itime));
     end
    
     if(sta_test==4)
     stas{1}(:,:,1,itime)=squeeze(fastClipSTA(:,:,itime));
     end
end
movieL=size(mov_modify_new,3);
nTrials=50;

movieTest = mov_orig;     %mov_orig,mov_modify_new
Response_test_movie_script
spksOrig= spksTrialsmall;

movieTest = mov_modify_new;%mov_orig,mov_modify_new
Response_test_movie_script
spksModi=spksTrialsmall;

figure;
subplot(2,1,1)
plotSpikeRaster(logical(spksOrig'),'PlotType','vertline');
title('Original Movie spike Raster');

subplot(2,1,2)
plotSpikeRaster(logical(spksModi'),'PlotType','vertline');
title('Null space movie spike Raster');
subplot(2,1,2)



