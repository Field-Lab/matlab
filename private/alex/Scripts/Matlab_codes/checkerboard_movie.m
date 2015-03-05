
cd('D:\Adapt_experiment\Stimuli_files')
clear all
% parameters
imageSize=600;
checkerSize=20;
time2pres=2;  %time in s
pres_cycles=time2pres*60;

mu=128;
sigma=20;
a=sigma*randn((imageSize/checkerSize)^2*pres_cycles,1)+mu;
a=round(a);
a(a<0)=0;
a(a>255)=255;
m=unique(a)';
col2take=reshape(a,pres_cycles,(imageSize/checkerSize)^2);
clear a

h=figure
for fr=1:pres_cycles
    c=zeros(imageSize,imageSize,3);
    cnt=1;
    for i=1:checkerSize:imageSize
        for j=1:checkerSize:imageSize
            c(i:i+checkerSize-1,j:j+checkerSize-1,1:3)=col2take(fr,cnt);
            cnt=cnt+1;
        end
    end
    c=uint8(c);
    image(c);
    F(fr)=getframe(h);
%     nn=['0000',int2str(fr)];
%     nn=nn(end-3:end);
%     imwrite(c,['/Users/alexth/Desktop/checkerboard_movie_mu128_sigma30/frame_',nn,'.png'],'png')
end

movie2avi(F, '/Users/alexth/Desktop/test.avi','compression','none','fps',60)

