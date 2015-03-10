function [STA,STC] = calculate_sta_stc(spksGen,mov,Filtlen)


Filtdim1=size(mov,1);
Filtdim2=size(mov,2);

binnedResponses = spksGen;
Len = length(spksGen);

% My own STA code 

STA=zeros(Filtdim1,Filtdim2,3,Filtlen);


for ibin=Filtlen+1:Len
   if(mod(ibin,1000)==1)
       ibin
   end
STA=STA+mov(:,:,:,ibin:-1:ibin-Filtlen+1)*binnedResponses(ibin);
end

STA=STA/sum(binnedResponses);

%% Process STA
STAs=sum(STA,3);
[sig_stixels, params, rf_strength_out] = significant_stixels(STAs,'select','thresh','thresh',2);
[r,c]=find(sig_stixels);
r=[min(r):max(r)];
c=[min(c):max(c)];
STA=squeeze(sum(STA,3));

 figure
 for itime=1:Filtlen
 imagesc(squeeze((STA(:,:,itime)))');colormap gray
 caxis([min(STA(:)),max(STA(:))]);
 colorbar
 pause(1/120)
 end

%% STC calculation 
display('STC calculation');
rSTA=STA(r,c,:); % Doubt!

 figure
 for itime=1:Filtlen
 imagesc(squeeze((rSTA(:,:,itime)))');colormap gray
 caxis([min(STA(:)),max(STA(:))]);
 colorbar
 pause(1/120)
 end

Filtdim1=length(r);
Filtdim2=length(c);

xxreSTA = reshape(rSTA,[Filtdim1*Filtdim2*Filtlen,1]);
xxreSTA=xxreSTA/norm(xxreSTA(:));


binnedResponsesTrial=binnedResponses;
framesValid = find(binnedResponsesTrial>0);

% flatten movie
mov=mov(r,c,:,:);
mov=squeeze(sum(mov,3));
reSTC=zeros(Filtdim1*Filtdim2*Filtlen,Filtdim1*Filtdim2*Filtlen);
for iframe=framesValid'
    iframe
  xx= reshape(mov(:,:,iframe:-1:iframe-Filtlen+1),[Filtdim1*Filtdim2*Filtlen,1]);
%  mask!

%  xx=  reshape(mov_new2(:,:,iframe:-1:iframe-Filtlen+1),[Filtdim1*Filtdim2*Filtlen,1]);
% Don't mask!

xxnew = (xx-((xx'*xxreSTA)*xxreSTA));
  
reSTC=reSTC+ xxnew*xxnew'*binnedResponsesTrial(iframe);
end


reSTC = reSTC / (sum(binnedResponses)-1);

[u,s,v]=svds(reSTC,100);
STC=reSTC;
figure;
plot(diag(s),'*')

figure;
for itime=1:30
subplot(2,2,1);
xu=reshape(u(:,1),[Filtdim1,Filtdim2,Filtlen]);
imagesc(xu(:,:,itime));
caxis([min(xu(:)),max(xu(:))]);
colormap gray

subplot(2,2,2);
xu=reshape(u(:,2),[Filtdim1,Filtdim2,Filtlen]);
imagesc(xu(:,:,itime));
caxis([min(xu(:)),max(xu(:))]);
colormap gray

subplot(2,2,3);
xu=reshape(u(:,3),[Filtdim1,Filtdim2,Filtlen]);
imagesc(xu(:,:,itime));
caxis([min(xu(:)),max(xu(:))]);
colormap gray

subplot(2,2,4);
xu=reshape(u(:,4),[Filtdim1,Filtdim2,Filtlen]);
imagesc(xu(:,:,itime));
caxis([min(xu(:)),max(xu(:))]);
colormap gray

itime
pause(1);
end


end