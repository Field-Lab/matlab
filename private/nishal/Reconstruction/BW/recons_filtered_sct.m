% Try to reconstruct the filtered cell activity
% STA not so good .. 
mov_norm=(mov-0.5);
filtered_act=zeros(finalLen,1);
sta_imp=squeeze(sta).*repmat(full(sig_stixels),[1,1,30]);
filter_len_use=5
for itime=30:finalLen
filtered_act(itime)=sum(sum(sum(mov_norm(:,:,itime-filter_len_use+1:itime).*sta_imp(:,:,30-filter_len_use+1:30))));
end

plot(filtered_act)

%%
delay=10;
negdelay=0;
trainsampleLen=round(finalLen*0.75);
noCells=nCells;

mov_recons = filtered_act-mean(filtered_act);
vecLen=(delay+negdelay)*noCells;
q=zeros(vecLen,1);
P=zeros(vecLen,vecLen);
nSamples=trainsampleLen-delay;
A=zeros(nSamples,vecLen);
b=zeros(nSamples,1);

for iLen=40:1:nSamples%1800*120
    iLen
    a=[];% add a constant
    for icell=1:noCells
    a=[a;double(spk_coll{icell}(iLen-negdelay:1:iLen+delay-1))];
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
    a=[];
    for icell=1:noCells
    a=[a;double(spk_coll{icell}(iLen-negdelay:1:iLen+delay-1))];
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