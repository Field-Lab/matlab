% Checkerboard flicker

cd('D:\Adapt_experiment\Stimuli_files')
clear all
% parameters
imageSize=600;
checkerSize=20;
time2pres=15;  %time in s
pres_cycles=time2pres*60;

mu=30;
sigma=9;
a=sigma*randn(900*pres_cycles,1)+mu;
a=round(a);
a(a<0)=0;
a(a>255)=255;
m=unique(a)';
col2take=reshape(a,pres_cycles,900);
high_contr=col2take;
clear a

corr_factor=0.2; % multiply sigma by this value

sigma_corr_factor=round(255*corr_factor);
mu_corr_factor=(mu-mu*corr_factor)/255;



% create templates textures
steps=imageSize/checkerSize;
[w, screenRect] = Screen(0,'OpenWindow',[],[100 100 1400 900],32); % open window, get screen size
scr=zeros(1,255);
% make texture of color 1 and 2 - once
for i=m
    imageArray=zeros(checkerSize,checkerSize)+i;
    imageArray=uint8(cat(3,imageArray,imageArray,imageArray));
    scr(i+1)=Screen('MakeTexture',w,imageArray);
end
% Screen('CloseAll');


%create coordinates map - once
center=screenRect(3:4)/2;
c=zeros(4,steps*steps);
cnt=1;
for i=1:steps
    for j=1:steps
        c(1,cnt)=center(1)+(i-steps/2-1)*checkerSize;
        c(2,cnt)=center(2)+(j-steps/2-1)*checkerSize;
        cnt=cnt+1;
    end
end
c(3,:)=c(1,:)+checkerSize;
c(4,:)=c(2,:)+checkerSize;
% run presentation
b=zeros(1,pres_cycles);
h=0;h1=0;
% a=round(rand(pres_cycles,steps*steps))+1;
% FlipTimestamp=zeros(1,pres_cycles);
FlipTimestamp=zeros(1,4200);

glsl1 = MakeTextureDrawShader(w, 'SeparateAlphaChannel',[mu_corr_factor*ones(1,3),0]);
tic
Screen('DrawTextures',w,scr(col2take(1,:)+1),[],c);
[h h1 FlipTimestamp(1)]=Screen('Flip',w,0,1,0);
tt=FlipTimestamp(1);
i=2;k=2;
for j=1:8
    if mod(j,2)==1
        while FlipTimestamp(i-1)-tt<5
            Screen('DrawTextures',w,scr(col2take(k,:)+1),[],c)
            [h h1 FlipTimestamp(i)]=Screen('Flip',w,0,1,0);
            i=i+1;
            k=k+1;
        end
        tt=FlipTimestamp(i-1);
        k=1;
    else
        while FlipTimestamp(i-1)-tt<5
            Screen('DrawTextures',w,scr(col2take(k,:)+1),[],c,[],[],[],sigma_corr_factor*ones(1,3),glsl1);
            [h h1 FlipTimestamp(i)]=Screen('Flip',w,0,1,0);
            i=i+1;
            k=k+1;
        end
        k=1;
        tt=FlipTimestamp(i-1);
    end
end
toc
Screen('CloseAll');



figure
% plot(diff(b(2:i-1)))
hold on
plot(diff(FlipTimestamp(2:i-1)),'r')




plot(diff(b(1:i-1)))
hold on
plot(diff(FlipTimestamp(1:i-1)),'r')
w
sum(diff(b(1:i-1)))

figure
dd=dir(['F:\20111103\flickers_timing\']);
for i=1:length(dd)-2
    load(['F:\20111103\flickers_timing\',dd(i+2).name])
    plot(diff(FlipTimeStamps))
    hold on
end



processed=double(processed(:,:,1));
mean(processed(:))
std(processed(:))

figure(2)
cnt=1;
clear bb gg
b=processed(:);
k=length(b);
for i=min(b):max(b)
    gg(cnt)=i;
    bb(cnt)=sum(b==i)/k;
    cnt=cnt+1;
end
% hold on
plot(gg,bb,'b')



not_processed=double(not_processed(:,:,1));
mean(not_processed(:))
std(not_processed(:))


cnt=1;
clear bb gg
b=not_processed(:);
k=length(b);
for i=min(b):max(b)
    gg(cnt)=i;
    bb(cnt)=sum(b==i)/k;
    cnt=cnt+1;
end
hold on
plot(gg,bb,'k')



% create full field sequences
mu=30;
sigma=9;
a=sigma*randn(900,1)+mu;
a=round(a);
a(a<0)=0;
a(a>255)=255;

mu=30;
sigma=1.8;
a=sigma*randn(900,1)+mu;
a=round(a);
a(a<0)=0;
a(a>255)=255;




