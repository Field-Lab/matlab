
%% Parameters
% Assume a grid of 320x640 


gridX=320;
gridY=640;

frameRate=120;
interval=2;
time=10; % in seconds
movieLen=time *frameRate/interval;

coarseinterval=10;
fineSz=4;

coarsemovieLen = time*frameRate/interval; 
CoarseSz = 16;

coarseType='checkerboard';
coarseframeRate = 1;

coarseScaleMax=0.24;
coarseScaleMin=0.12;

fineScale=0.12;
%% coarseMesh
coarseMesh=zeros(gridX/CoarseSz,gridY/CoarseSz);
%
if(strcmp(coarseType,'random'))
coarseMesh = (rand(size(coarseMesh))>0.5 );
end

if(strcmp(coarseType,'checkerboard'))
icnt=0;
for icntx=1:size(coarseMesh,1)
    for icnty=1:size(coarseMesh,2)
       icnt=icnt+1;
       if(rem(icnt+icntx,2)==1)
    coarseMesh(icntx,icnty)=1;
       end
    end
end
end

MovieCoarse=zeros(gridX/CoarseSz,gridY/CoarseSz,coarsemovieLen);
MovieCoarse = MovieCoarse+2*((rand(size(MovieCoarse))>0.5)-0.5).*repmat(coarseScaleMin + (coarseScaleMax - coarseScaleMin)*coarseMesh,[1,1,movieLen]);


CoarseMovieFine = zeros(gridX/fineSz,gridY/fineSz,movieLen);
for icntx=1:size(CoarseMovieFine,1)
    for icnty=1:size(CoarseMovieFine,2)
        for itime=1:movieLen
    CoarseMovieFine(icntx,icnty,itime)=MovieCoarse(floor((icntx-1)*fineSz/CoarseSz)+1,floor((icnty-1)*fineSz/CoarseSz)+1,floor((itime-1)*interval/coarseinterval)+1);
        end
    end
end

%% coarseMovie

MovieFine=fineScale*2*((rand(gridX/fineSz,gridY/fineSz,movieLen)>0.5)-0.5);
Movie_full = MovieFine+CoarseMovieFine;

% 
% writerObj = VideoWriter('/Volumes/Lab/Users/bhaishahster/perturb.avi','Grayscale AVI');
% writerObj.FrameRate=120;
% open(writerObj);
% writeVideo(writerObj,Movie_full+0.5);
% close(writerObj);


figure;
for itime=1:movieLen
    subplot(1,2,1);
    imagesc(CoarseMovieFine(:,:,itime)');
    colormap gray
    axis image
    
    subplot(1,2,2);
    imagesc(Movie_full(:,:,itime)');
    colormap gray
    axis image
    pause(1)
end