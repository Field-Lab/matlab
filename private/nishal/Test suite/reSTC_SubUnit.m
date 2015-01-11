
step=30;
nTrials=1;
% My own STA code 
reSTA=zeros(Filtdim1,Filtdim2,Filtlen);
icnt=0;
for itrial=1:nTrials
for iframe=59:step:movie_new_len
    if(binnedResponses(iframe,itrial)>0)
reSTA=reSTA+repmat(mask,[1,1,30]).*mov_new2(:,:,iframe:-1:iframe-Filtlen+1)*binnedResponses(iframe,itrial);
icnt=icnt+binnedResponses(iframe,itrial);
    end
end
end
reSTA=reSTA/icnt;


figure
 for itime=1:Filtlen
 imagesc(squeeze(reSTA(:,:,itime)));colormap gray
 caxis([min(reSTA(:)),max(reSTA(:))]);
 colorbar
 pause(1/120)
 end
 
 % Relation between STA and re-STA

 % Histograms
 figure
 [X,C]=hist(STA(:),20);
 [X2,C2]=hist(reSTA(:),20);
 plotyy(C,X,C2,X2);
 title('Histograms of STA and re-STA')
 legend('STA','reSTA')
 
 % Correlation between STA and re-STA
 corrSTA=STA.*reSTA;
 corr_STAr=sum(corrSTA(:))/(norm(STA(:))*norm(reSTA(:))) 
 disp('Obviously, Correlation between STA and reSTA too low as its in null space.. ')
 
 %% Re- STC
 
% My own STA code 



reSTC=zeros(Filtdim1*Filtdim2*Filtlen);
indx=[59:1:movie_new_len];
for iTrial=1:nTrials

binnedResponsesTrial=binnedResponses(:,iTrial);



framesValid = indx(binnedResponsesTrial(indx)>0);

for iframe=framesValid
    iframe
  xx=  reshape(repmat(mask,[1,1,30]).*mov_new2(:,:,iframe:-1:iframe-Filtlen+1),[Filtdim1*Filtdim2*Filtlen,1]);
  
reSTC=reSTC+ xx*xx'*binnedResponsesTrial(iframe);
end
end
yy=reshape(reSTA,[Filtdim1*Filtdim2*Filtlen,1]);
reSTC=reSTC / (sum(binnedResponses(:))-1) - yy*yy'* (sum(binnedResponses(:)))/(sum(binnedResponses(:))-1);

[u,s,v]=svds(reSTC,100);
figure;
plot(diag(s),'*')
%%
 
figure;
for isubunit=1:nSubunits
subplot(2,2,isubunit);
imagesc(subunits{isubunit}(:,:,4));
colormap gray
colorbar
end

uSq=cell(4,1); 
for isub=1:4
uSq{isub}=reshape(u(:,isub),[Filtdim1,Filtdim2,Filtlen])%.*repmat(mask,[1,1,Filtlen]);
end


 figure
 
 for itime=4%1:Filtlen
 for isub=1:4
     subplot(2,2,isub);
     imagesc(squeeze(uSq{isub}(:,:,itime)));
 colormap gray
 caxis([min(uSq{isub}(:)),max(uSq{isub}(:))])
 colorbar
 %pause
 end
 end
 
%  figure
%  for itime=1:Filtlen
%      itime
%      for iidx=1:36
%      subplot(6,6,iidx)
%      uSq=reshape(u(:,iidx),[Filtdim1,Filtdim2,Filtlen]).*repmat(mask,[1,1,Filtlen]);
%      imagesc(squeeze(uSq(:,:,itime)));
%      colormap gray
%      %caxis([min(u(:)),max(u(:))]);
%       caxis([min(uSq(:)),max(uSq(:))]);
%      colorbar
%      end
% pause
%  end
%  
 
%uSq=reshape(u(:,72),[Filtdim1,Filtdim2,Filtlen]).*repmat(mask,[1,1,Filtlen]);
% cell_resp=zeros(1,nSubunits);
% 
% for isubunit=1:nSubunits
%    
% cell_resp(:,isubunit)=sum(uSq(:).*squeeze(subunits{isubunit}(:)));
% end
% cell_resp
% % 
% totalInput=0*cell_resp(:,1);
% for isubunit=1:nSubunits
% totalInput=totalInput+subunitWeights(isubunit)*f(cell_resp(:,isubunit));
% end
% 
% totalOutput=N(totalInput)
