function response = calculate_sta_ts(mov_params,response,sta_params,cell_params)

spksGen=response.spksGen;
binsLen = length(spksGen);
bin_to_frame=response.mov_frame_number;
mov=mov_params.mov;
Filtlen=sta_params.Filtlen;
useTrial= sta_params.useTrial;
binsPerFrame= cell_params.binsPerFrame;
Filtdim1=size(mov,1);
Filtdim2=size(mov,2);

binnedResponses = spksGen;
% My own STA code 
STA=zeros(Filtdim1,Filtdim2,3,Filtlen);
binnedFramesResponse=zeros(max(bin_to_frame),1);
ibin=1;
for iframe=1:max(bin_to_frame)
    
binnedFramesResponse(iframe)=sum(binnedResponses(useTrial,ibin:ibin+binsPerFrame-1));

ibin=ibin+binsPerFrame;
end

idx=1:length(binnedFramesResponse);
frames = idx(binnedFramesResponse==1);
for iframe=frames(frames>Filtlen+1)
   if(mod(iframe,1000)==1)
       iframe
   end
STA=STA+mov(:,:,:,iframe:-1:iframe-Filtlen+1)*binnedFramesResponse(iframe);
end

STA=STA/sum(binnedResponses(useTrial,:));


STA=squeeze(sum(STA,3));

 figure
 for itime=1:Filtlen
 imagesc(squeeze((STA(:,:,itime)))');colormap gray
 caxis([min(STA(:)),max(STA(:))]);
 axis image
 colorbar
 pause(1/120)
 end
 response.analyse.STA=STA;
 
end