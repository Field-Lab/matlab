clear

pathway='/Users/alexth/Desktop/old_stuff/20140121b/';
words=['#0009';'#0010';'#0072';'#0102';'#0164';'#0165';'#0166';'#0228';'#0229';'#0230';'#0292';'#0293';'#0294'];

for cnt=9:size(words,1)
    sta=0;
    codeWord=words(cnt,:);
    tic
    sta=sta_from_CB(codeWord, pathway, sta);
    toc
    save([pathway,'sta_',codeWord,'.mat'],'sta');
end



load([pathway,'sta_#0228'])

for unitN=5
    figure
    sta_tmp=sta(:,:,unitN)';
    sta_tmp=reshape(sta_tmp,40,40,500);
    cnt=1;
    for i=1:20:500
        subplot(5,5,cnt)
        colormap gray
        imagesc(sta_tmp(:,:,i))
        cnt=cnt+1;
    end
end




load([pathway,'sta_#0228'])
tmp=sta;
load([pathway,'sta_#0229'])
tmp=tmp+sta;
load([pathway,'sta_#0230'])
tmp=tmp+sta;

for unitN=5
    figure
    sta_tmp=tmp(:,:,unitN)';
    sta_tmp=reshape(sta_tmp,40,40,500);
    cnt=1;
    for i=1:20:500
        subplot(5,5,cnt)
        colormap gray
        imagesc(sta_tmp(:,:,i))
        cnt=cnt+1;
    end
end



load([pathway,'sta_#0072'])
tmp1=sta;
load([pathway,'sta_#0102'])
tmp1=tmp1+sta;



load([pathway,'sta_#0164'])
tmp2=sta;
load([pathway,'sta_#0165'])
tmp2=tmp2+sta;
load([pathway,'sta_#0166'])
tmp2=tmp2+sta;



for unitN=5
    figure
    a=tmp(:,:,unitN)'-min(min(tmp(:,:,unitN)));
    a=a/max(a(:));
    b=tmp1(:,:,unitN)'-min(min(tmp1(:,:,unitN)));
    b=b/max(b(:));
    c=tmp2(:,:,unitN)'-min(min(tmp2(:,:,unitN)));
    c=c/max(c(:));
    
    sta=zeros(40,40,500,3);
    sta(:,:,:,1)=reshape(a,40,40,500);
    sta(:,:,:,2)=reshape(b,40,40,500);
    sta(:,:,:,3)=reshape(c,40,40,500);
    sta=permute(sta,[1 2 4 3]);
    
    sta=zeros(40,40,500,3);
    sta(:,:,:,1)=0;
    sta(:,:,:,2)=reshape(a,40,40,500);
    sta(:,:,:,3)=0;
    sta=permute(sta,[1 2 4 3]);

        cnt=1;
    for i=125:25:500
        subplot(4,4,cnt)
        imagesc(sta(:,:,:,i))
        cnt=cnt+1;
    end
end




