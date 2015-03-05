
cd('D:\Adapt_experiment\Stimuli_files')
clear all
% parameters
imageSize=600;
checkerSize=20;
color1=20;
color2=200;
time2pres=10;  %time in s

pres_cycles=time2pres*60;
for i=1:20
    checkers=logical(round(rand(pres_cycles,steps*steps)));
    save(['\\134.2.181.58\ag_muench\user\alexandra\checkers\checkers',int2str(i)],'checkers');
end


% create templates textures
steps=imageSize/checkerSize;
imageArray=zeros(checkerSize,checkerSize)+color1;
imageArray2=zeros(checkerSize,checkerSize)+color2;
imageArray=uint8(cat(3,imageArray,imageArray,imageArray));
imageArray2=uint8(cat(3,imageArray2,imageArray2,imageArray2));
[w, screenRect] = Screen(0,'OpenWindow',[],[],32); % open window, get screen size

% make texture of color 1 and 2 - once
scr(1)=Screen('MakeTexture',w,imageArray);
scr(2)=Screen('MakeTexture',w,imageArray2);
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
a=round(rand(pres_cycles,steps*steps))+1;
FlipTimestamp=zeros(1,pres_cycles);
FlipTimestamp(1)=1;

Screen('DrawTextures',w,scr(a(1,:)),[],c);
%     Screen('WaitBlanking', w,1);
[h h1 FlipTimestamp(1)]=Screen('Flip',w,0,1,0);
tic
i=2;
while FlipTimestamp(i-1)-FlipTimestamp(1)<time2pres-2/60
    Screen('DrawTextures',w,scr(a(i,:)),[],c);
%     Screen('WaitBlanking', w,1);
    [h h1 FlipTimestamp(i)]=Screen('Flip',w,0,1,0);
    b(i)=toc;
    i=i+1;
end
toc
Screen('CloseAll');
FlipTimestamp=FlipTimestamp(1:i-1);
b=b(2:i-1);
plot(FlipTimestamp-FlipTimestamp(1),'Marker','x')
plot(diff(FlipTimestamp-FlipTimestamp(1)),'Marker','x')
plot(diff(b));
% close all
% plot(b)
% hold on
% plot(d,'color','r')
tic
load('C:\Documents and Settings\atikidzhi\Desktop\b.mat')
toc

