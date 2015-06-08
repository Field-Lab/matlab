function [WNSTA,WNSTC,WN_uSq,h] = compute_STA_STC(binnedResponses,mov2)
Filtdim1 = size(mov2,1);
movieLen = size(mov2,2);
mov=zeros(Filtdim1,1,movieLen);
mov(:,1,:)=mov2;
Filtdim2 = 1;

Filtlen=1;
% My own STA code 
STA=zeros(Filtdim1,Filtdim2,Filtlen);
useTrial=1;
for iframe=30:movieLen
STA=STA+mov(:,:,iframe:-1:iframe-Filtlen+1)*binnedResponses(iframe(:,useTrial));
end
STA=STA/sum(binnedResponses(:,useTrial));

% 
% figure
%  for itime=[1:Filtlen]
%  imagesc(squeeze(STA(:,:,itime)));colormap gray
%  caxis([min(STA(:)),max(STA(:))]);
%  colorbar
%  pause(1/120)
%  end
 
movie_new_len = size(mov,3);
mov_new2=mov;
reSTC_SubUnit_subtractSTA % Calculate STA and STC .. 
%reSTC_SubUnit; % Do not subtract/remove STA
WNSTA=reSTA;
WNSTC=reSTC;

WN_uSq=uSq;


end