
%finalLen from before
%noCells from before

delay=3;
trainsampleLen=round(finalLen*0.75);
noCells=nCells;

filter_collection=cell(length(sig_pix),1);

for sig_pix_idx=1:length(sig_pix);
pixX=sig_pix(sig_pix_idx,1);
pixY=sig_pix(sig_pix_idx,2);

pixCol=1:1;
pixX
pixY
pixCol

mov_recons = squeeze(mov(pixX,pixY,pixCol,:))-0.5;
vecLen=delay*noCells+1;
q=zeros(vecLen,1);
P=zeros(vecLen,vecLen);
nSamples=trainsampleLen-delay;
A=zeros(nSamples,vecLen);
b=zeros(nSamples,1);

for iLen=1:1:nSamples%1800*120
    iLen
    a=[1];% add a constant
    for icell=1:noCells
    a=[a;double(spk_coll{icell}(iLen:1:iLen+delay-1))];
    end
    a=full(a);
    y=mov_recons(iLen);

    A(iLen,:)=a';
    b(iLen)=y;
    
    q=q+y*a;
    P=P+a*a';

end

 
filter=P\q;   
%filter2=(A\(A'\(A'*b)));

filter_use=filter;

% Make prediction
mov_pred=0*mov_recons;
recons_idx=1:trainsampleLen%trainsampleLen+1:finalLen-delay%1:trainsampleLen%trainsampleLen+1:finalLen;%1:trainsampleLen%
for iLen=recons_idx %1800*120
    iLen
    a=[1];
    for icell=1:noCells
    a=[a;double(spk_coll{icell}(iLen:1:iLen+delay-1))];
    end
    a=full(a);
    
   mov_pred(iLen)= filter_use'*a;
    
end

% % m1= (mov_recons(recons_idx(1)+100:iLen));
% % m2= (mov_pred(recons_idx(1)+100:iLen));
% % figure;
% % stairs(m1,'b');
% % hold on
% % stairs(m2,'r');
% % %hold on;
% % ylim([-1,1]);
%stairs(double(mov_pred(recons_idx(1)+100:recons_idx(1)+200)>0)-0.5,'g');
%xlim([1,100])
% correlation 

% figure;
% %filter_bank{pixX,pixY,pixCol}.filter=filter_use;
% %filter_bank{pixX,pixY,pixCol}.correlation=corr(mov_pred(recons_idx),mov_recons(binFrames(recons_idx)));
% corr(mov_pred(recons_idx(1)+100:recons_idx(1)+2000),mov_recons(recons_idx(1)+100:recons_idx(1)+2000))
% corr(m1,m2)
% [r,lags]=xcorr(m1,m2);
% plot(lags,r);
%% See filters
%figure
  filter_len=length(q);
cellFilter_len=(filter_len-1)/nCells;
filter_coll=zeros(cellFilter_len,nCells);
z=filter_use;
for iF=1:noCells
    z_cut=z((iF-1)*cellFilter_len+2:iF*cellFilter_len+1);
filter_coll(:,iF)=z_cut;
end
%plot(filter)

filter_collection{sig_pix_idx}=filter_use;
end


%% 


mov_pred_full=zeros(size(mov,1),size(mov,2),length(mov_recons));
mov_recons_full=zeros(size(mov,1),size(mov,2),length(mov_recons));

for sig_pix_idx=1:length(sig_pix);
pixX=sig_pix(sig_pix_idx,1);
pixY=sig_pix(sig_pix_idx,2);

pixCol=1:1;
pixX
pixY
pixCol

filter_use=filter_collection{sig_pix_idx};

% Make prediction

recons_idx=trainsampleLen+1:finalLen-delay%1:trainsampleLen%trainsampleLen+1:finalLen-delay%1:trainsampleLen%trainsampleLen+1:finalLen;%1:trainsampleLen%
for iLen=recons_idx %1800*120
    iLen
    a=[1];
    for icell=1:noCells
    a=[a;double(spk_coll{icell}(iLen:1:iLen+delay-1))];
    end
    a=full(a);
    
   mov_pred(iLen)= filter_use'*a;
    
end

mov_recons = squeeze(mov(pixX,pixY,pixCol,:))-0.5;
m1= (mov_recons(recons_idx(1)+100:recons_idx(1)+500));
m2= (mov_pred(recons_idx(1)+100:recons_idx(1)+500));
figure;
stairs(m1,'b');
hold on
stairs(m2,'r');
%hold on;
ylim([-1,1]);
title(sprintf('sig pix idx %d',sig_pix_idx));
%stairs(double(mov_pred(recons_idx(1)+100:recons_idx(1)+200)>0)-0.5,'g');
%xlim([1,100])
% correlation 


mov_pred_full(pixX,pixY,:)=mov_pred;
mov_recons_full(pixX,pixY,:)=mov_recons;
end

%% 
figure
for idx=recons_idx(100:end-100)
subplot(2,1,1);
imagesc(mov_recons_full(:,:,idx));
colormap gray
caxis([-0.5,0.5]);
colorbar
title('Spatial')
subplot(2,1,2);
imagesc(mov_pred_full(:,:,idx));
colormap gray
caxis([-0.5,0.5]);
colorbar

pause

end

%% DOUBTFUL
wellFittedStixIdx = 1:16;%[2,3,6,7,10,11,12,14,15,16];
v=zeros(size(mov_recons_full,3),1);
u=zeros(size(mov_recons_full,3),1);
for iLen = 1:length(wellFittedStixIdx)
v=v+squeeze(mov_recons_full(sig_pix(wellFittedStixIdx(iLen),1),sig_pix(wellFittedStixIdx(iLen),2),:));
u=u+squeeze(mov_pred_full(sig_pix(wellFittedStixIdx(iLen),1),sig_pix(wellFittedStixIdx(iLen),2),:));
end

figure;
stairs(v,'b');
hold on
stairs(u,'r')
legend('Original','Reconstructed');


