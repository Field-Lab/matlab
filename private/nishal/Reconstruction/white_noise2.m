
%finalLen from before
%noCells from before

delay=20;
trainsampleLen=round(finalLen*0.75);
noCells=nCells;

pixX=15;
pixY=15;

pixCol=1:1;
pixX
pixY
pixCol

mov_recons = squeeze(mov(pixX,pixY,pixCol,:))-0.5;
q=zeros(delay*noCells+1,1);
P=zeros(delay*noCells+1,delay*noCells+1);
nSamples=trainsampleLen-delay;
A=zeros(nSamples,delay*noCells+1);
b=zeros(nSamples,1);

for iLen=1:nSamples%1800*120
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
recons_idx=trainsampleLen+1:finalLen-delay%1:trainsampleLen%trainsampleLen+1:finalLen;%1:trainsampleLen%
for iLen=recons_idx %1800*120
    iLen
    a=[1];
    for icell=1:noCells
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
%stairs(double(mov_pred(recons_idx(1)+100:recons_idx(1)+200)>0)-0.5,'g');
%xlim([1,100])
% correlation 

%filter_bank{pixX,pixY,pixCol}.filter=filter_use;
%filter_bank{pixX,pixY,pixCol}.correlation=corr(mov_pred(recons_idx),mov_recons(binFrames(recons_idx)));
corr(mov_pred(recons_idx(1)+100:recons_idx(1)+2000),mov_recons(recons_idx(1)+100:recons_idx(1)+2000))
corr(m1,m2)
[r,lags]=xcorr(m1,m2);
plot(lags,r);
%% See filters

  filter_len=length(q);
cellFilter_len=(filter_len-1)/nCells;
filter_coll=zeros(cellFilter_len,nCells);
z=filter_use;
for iF=1:noCells
    z_cut=z((iF-1)*cellFilter_len+2:iF*cellFilter_len+1);
filter_coll(:,iF)=z_cut;
end