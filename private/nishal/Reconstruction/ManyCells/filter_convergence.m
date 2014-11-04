
%finalLen from before
%noCells from before

delay=3;
noCells=nCells;
percent_list=[0.1,0.3,0.5,0.7,1,2,3,4,5,10,20,30,50,60,75,80];
for percent=percent_list;
    fraction=percent/100
trainsampleLen=round(finalLen*fraction);


filter_collection=struct('filter_use',[]);

mov_pred_full=zeros(size(mov,1),size(mov,2),size(mov,4));
mov_recons_full=zeros(size(mov,1),size(mov,2),size(mov,4));

for pixX=12:17;
    for pixY=12:17;

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
vecLen=delay*sig_cells;
q=zeros(vecLen,1);
P=zeros(vecLen,vecLen);
nSamples=trainsampleLen-delay;
A=zeros(nSamples,vecLen);
b=zeros(nSamples,1);

for iLen=1:1:nSamples%1800*120
    iLen;
    a=[];% add a constant
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
recons_idx=floor(finalLen*0.8):finalLen-delay;%1:trainsampleLen%trainsampleLen+1:finalLen;%1:trainsampleLen%
for iLen=recons_idx %1800*120
    iLen;
    a=[];
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


mov_pred_full(pixX,pixY,:)=mov_pred;
mov_recons_full(pixX,pixY,:)=mov_recons;

close all


filter_collection(pixX,pixY).filter_use=filter_use;
filter_collection(pixX,pixY).cell_list=cell_list;

    end
end
if(percent<1)
save(sprintf('/Volumes/Analysis/nishal/Reconstruction/recons_WN_On_percent_%f.mat',percent),'filter_collection','sta_data','mov_pred_full','mov_recons_full','percent')

else
    
save(sprintf('/Volumes/Analysis/nishal/Reconstruction/recons_WN_On_percent_%d.mat',percent),'filter_collection','sta_data','mov_pred_full','mov_recons_full','percent')
end
end



%%
recons_idx=[finalLen*0.8:finalLen-delay]; % Doubt 

h=fspecial('gaussian',4,2);

corr_data=[];
sq_err=[];
icnt=0;
for  percent=percent_list;
    icnt=icnt+1
    if(percent>=1)
    Data=load(sprintf('/Volumes/Analysis/nishal/Reconstruction/recons_WN_On_percent_%d.mat',percent'));
    else
       Data=load(sprintf('/Volumes/Analysis/nishal/Reconstruction/recons_WN_On_percent_%f.mat',percent'));
  
    end
    Data.mov_recons_full_smooth=zeros(size(Data.mov_recons_full));
for idx=floor(0.8*finalLen):finalLen-delay      
    Data.mov_recons_full_smooth(:,:,idx)=imfilter(Data.mov_recons_full(:,:,idx),h,'same').*(Data.mov_recons_full(:,:,idx)~=0);
end

corr_data(icnt)=corr(Data.mov_recons_full_smooth(:),Data.mov_pred_full(:));
sq_err(icnt)=sum(sum(sum((Data.mov_recons_full_smooth(:)-Data.mov_pred_full(:)).^2)));
end

figure;
plot(percent_list*duration*5*8.33/(100*1000*60),corr_data,'LineWidth',2);
title('Movie : BW-20-5-0.48-11111-30x30-60.35, 21741 frames');
xlabel('minutes')
ylabel('Test accuracy (Correlation)');

figure;
semilogy(percent_list*duration*5*8.33/(100*1000*60),sq_err,'LineWidth',2);
title('Movie : BW-20-5-0.48-11111-30x30-60.35, 21741 frames');
xlabel('minutes')
ylabel('Test accuracy (Sq. Error)');
