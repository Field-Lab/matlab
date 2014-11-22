function response = calculate_sta_ts(mov_params,response,sta_params)

spksGen=response.spksGen;
binsLen = length(spksGen);
bin_to_frame=response.mov_frame_number;
mov=mov_params.mov;
Filtlen=sta_params.Filtlen;

Filtdim1=size(mov,1);
Filtdim2=size(mov,2);

binnedResponses = spksGen;
% My own STA code 
STA=zeros(Filtdim1,Filtdim2,3,Filtlen);

for iframe=30:max(bin_to_frame)
    iframe
STA=STA+mov(:,:,:,iframe:-1:iframe-Filtlen+1)*sum(binnedResponses(bin_to_frame==iframe));
end

STA=STA/sum(binnedResponses);


STA=squeeze(sum(STA,3));

 figure
 for itime=1:Filtlen
 imagesc(squeeze((STA(:,:,itime)))');colormap gray
 caxis([min(STA(:)),max(STA(:))]);
 colorbar
 pause
 end
 response.analyse.STA=STA;
 
end