
%finalLen from before
%noCells from before

delay=3;
trainsampleLen=round(finalLen*0.75);
noCells=nCells;

filter_collection=struct('filter_use',[]);


for pixX=1:30;
    for pixY=1:30;

pixCol=1:1;
pixX
pixY
pixCol

close all
%% Find Significant cells
cell_list=[];
for icell=1:noCells
if(sta_data(icell).sig_stixels(pixX,pixY)==1)
   cell_list=[cell_list;icell]; 
end
end

sig_cells=length(cell_list)
cell_list


filter_collection(pixX,pixY).sig_cells=sig_cells;
if(sig_cells==0)
     continue
end
%%



mov_recons = squeeze(mov(pixX,pixY,pixCol,:))-0.5;
vecLen=delay*sig_cells+1;
q=zeros(vecLen,1);
P=zeros(vecLen,vecLen);
nSamples=trainsampleLen-delay;
A=zeros(nSamples,vecLen);
b=zeros(nSamples,1);

for iLen=1:1:nSamples%1800*120
    iLen;
    a=[1];% add a constant
    for icell=cell_list'
        icell;
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
recons_idx=trainsampleLen+1:finalLen-delay;%1:trainsampleLen%trainsampleLen+1:finalLen;%1:trainsampleLen%
for iLen=recons_idx %1800*120
    iLen;
    a=[1];
    for icell=cell_list'
    a=[a;double(spk_coll{icell}(iLen:1:iLen+delay-1))];
    end
    a=full(a);
    
   mov_pred(iLen)= filter_use'*a;
    
end

m1= (mov_recons(recons_idx(1)+100:iLen));
m2= (mov_pred(recons_idx(1)+100:iLen));
figure;
stairs(m1,'b');
hold on
stairs(m2,'r');
%hold on;
ylim([-1,1]);
title(sprintf('Pix %d %d',pixX,pixY));
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
figure
  filter_len=length(q);
cellFilter_len=(filter_len-1)/sig_cells;
filter_coll=zeros(cellFilter_len,sig_cells);
z=filter_use;
for iF=1:sig_cells
    z_cut=z((iF-1)*cellFilter_len+2:iF*cellFilter_len+1);
filter_coll(:,iF)=z_cut;
end
plot(filter_coll)
title(sprintf('Pix %d %d',pixX,pixY));

filter_collection(pixX,pixY).filter_use=filter_use;
filter_collection(pixX,pixY).cell_list=cell_list;

    end
end


%% 


mov_pred_full=zeros(size(mov,1),size(mov,2),length(mov_recons));
mov_recons_full=zeros(size(mov,1),size(mov,2),length(mov_recons));



for pixX=1:30;
    for pixY=1:30;
close all
        if(filter_collection(pixX,pixY).sig_cells==0)
        continue;
        end
        
pixCol=1:1;
pixX
pixY
pixCol


filter_use=filter_collection(pixX,pixY).filter_use;
cell_list=filter_collection(pixX,pixY).cell_list;
% Make prediction

%1:trainsampleLen%trainsampleLen+1:finalLen-delay%1:trainsampleLen%trainsampleLen+1:finalLen;%1:trainsampleLen%
recons_idx=trainsampleLen+1:finalLen-delay;%1:trainsampleLen%trainsampleLen+1:finalLen;%1:trainsampleLen%
for iLen=recons_idx %1800*120
    iLen;
    a=[1];
    for icell=cell_list'
    a=[a;double(spk_coll{icell}(iLen:1:iLen+delay-1))];
    end
    a=full(a);
    
   mov_pred(iLen)= filter_use'*a;
    
end

mov_recons = squeeze(mov(pixX,pixY,pixCol,:))-0.5;
m1= (mov_recons(recons_idx(1)+100:recons_idx(end)));
m2= (mov_pred(recons_idx(1)+100:recons_idx(end)));
figure;
stairs(m1,'b');
hold on
stairs(m2,'r');
%hold on;
ylim([-1,1]);
title(sprintf('Pix %d %d',pixX,pixY));
%stairs(double(mov_pred(recons_idx(1)+100:recons_idx(1)+200)>0)-0.5,'g');
%xlim([1,100])
% correlation 


mov_pred_full(pixX,pixY,:)=mov_pred;
mov_recons_full(pixX,pixY,:)=mov_recons;
end
end

%% 
figure
for idx=recons_idx(100:end-100)
subplot(3,2,1);
imagesc(mov_recons_full(:,:,idx));
colormap gray
axis image
caxis([-0.5,0.5]);
colorbar
%title('Spatial')

subplot(3,2,2);
imagesc(mov_pred_full(:,:,idx));
colormap gray
caxis([-0.5,0.5]);
colorbar
axis image

image_filter=[1/9,1/9,1/9;
             1/9,1/9,1/9;
              1/9,1/9,1/9];

subplot(3,2,3);
ppxx=imfilter(mov_recons_full(:,:,idx),image_filter,'same').*(mov_recons_full(:,:,idx)~=0);
imagesc(ppxx);
colormap gray
axis image
caxis([-0.5,0.5]);
colorbar

subplot(3,2,4);
imagesc(0.5*double(mov_pred_full(:,:,idx)>0) - 0.5*double(mov_pred_full(:,:,idx)<-0) );
colormap gray
caxis([-0.5,0.5]);
colorbar
axis image

image_filter=[1/9,1/9,1/9;
             1/9,1/9,1/9;
              1/9,1/9,1/9];

subplot(3,2,6);
ppxx=imfilter(mov_recons_full(:,:,idx),image_filter,'same').*(mov_recons_full(:,:,idx)~=0);
imagesc(0.5*double(ppxx>0) - 0.5*double(ppxx<-0));
colormap gray
axis image
caxis([-0.5,0.5]);
colorbar

pause

end

%% DOUBTFUL
% wellFittedStixIdx = 1:16;%[2,3,6,7,10,11,12,14,15,16];
% v=zeros(size(mov_recons_full,3),1);
% u=zeros(size(mov_recons_full,3),1);
% for iLen = 1:length(wellFittedStixIdx)
% v=v+squeeze(mov_recons_full(sig_pix(wellFittedStixIdx(iLen),1),sig_pix(wellFittedStixIdx(iLen),2),:));
% u=u+squeeze(mov_pred_full(sig_pix(wellFittedStixIdx(iLen),1),sig_pix(wellFittedStixIdx(iLen),2),:));
% end
% 
% figure;
% stairs(v,'b');
% hold on
% stairs(u,'r')
% legend('Original','Reconstructed');
% 
% % High correlation between total input and reconstruction .. 
% corr(u,v)
% 
% figure;
% plot(xcorr(u,v))