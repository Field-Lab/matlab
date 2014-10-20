
%%
% Fit STA
addpath(genpath('~/Dropbox/Lab/Development/matlab-standard/code')); 
addpath(genpath('~/Dropbox/Lab/Development/matlab-standard/private/nishal/fwdfittingfunctions'));
fit_info = fit_sta(sta_data002)
%datarun = compute_sta_fits(datarun, 'all');

%params=datarun.matlab.sta_fits{matlab_id}.initial_params;
params=fit_info.initial_params;
full_fit = sta_fit_function(params);

ssta=sta_data002;

newSTAFitSz=30;
sta_rsz=zeros(newSTAFitSz,newSTAFitSz,3,30);
for itime=1:30
x=imresize(full_fit(:,:,1,itime),[newSTAFitSz,newSTAFitSz]);
sta_rsz(:,:,1,itime)=x(end:-1:1,end:-1:1);
x=imresize(full_fit(:,:,2,itime),[newSTAFitSz,newSTAFitSz]);
sta_rsz(:,:,2,itime)=x(end:-1:1,end:-1:1);
x=imresize(full_fit(:,:,3,itime),[newSTAFitSz,newSTAFitSz]);
sta_rsz(:,:,3,itime)=x(end:-1:1,end:-1:1);

end

figure;
for itime=30:-1:1
    itime
    subplot(3,1,1);
    imagesc(squeeze(sum(ssta(:,:,:,itime),3)));
    colormap gray
    axis image
    colorbar
    caxis([min(ssta(:)), max(ssta(:))]);

    
    subplot(3,1,2);
    imagesc(squeeze(sum(full_fit(:,:,:,itime),3)));
    colormap gray
    axis image
    colorbar 
    caxis([min(ssta(:)), max(ssta(:))]);
% 
    subplot(3,1,3);
    imagesc(squeeze(sum(sta_rsz(:,:,:,itime),3)));
    colormap gray
    axis image
    colorbar 
    caxis([min(ssta(:)), max(ssta(:))]);
  
    pause
    
end


%%
% Try to reconstruct the filtered cell activity
% STA not so good .. 

% WHAT TO DO !!?? MOVIE AT COARSER RESOLUTION ..
mov_norm=(mov-0.5);
filtered_act=zeros(finalLen,1);
sta_imp=squeeze(sum(sta_rsz(:,:,:,:),3));
% 
sta_bin_sz=100/(distinctImagesbetweenTriggers*subdivideRefresh);
sta_binned=zeros(30,30,30/sta_bin_sz);
icnt=0;
for itime=1:sta_bin_sz:30
    itime
    icnt=icnt+1;
    sta_binned(:,:,icnt)=sum(sta_imp(:,:,itime:itime+sta_bin_sz-1),3);
end

filter_len_use=2

for itime=30:finalLen
filtered_act(itime)=sum(sum(sum(mov_norm(:,:,itime-filter_len_use+1:itime).*sta_binned(:,:,filter_len_use:-1:1))));
end

filtered_act=filtered_act*10;
stairs(filtered_act)

%%
delay=3;
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
%ylim([-1,1]);