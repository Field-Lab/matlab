
%% Parameters
% Assume a grid of 320x640 


gridX=320;
gridY=640;

frameRate=120;
interval=2;
movieLen=10 *frameRate/interval;

CoarseSz=16;
fineSz=4;

coarseType='checkerboard';

coarseScaleMax=0.24;
coarseScaleMin=0.06;

fineScale=0.24;
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

MovieCoarse=zeros(gridX/CoarseSz,gridY/CoarseSz,movieLen);
MovieCoarse = MovieCoarse+2*((rand(size(MovieCoarse))>0.5)-0.5).*repmat(coarseScaleMin + (coarseScaleMax - coarseScaleMin)*coarseMesh,[1,1,movieLen]);

%% coarseMovie

MovieFine=fineScale*2*((rand(gridX/fineSz,gridY/fineSz,movieLen)>0.5)-0.5);
Movie_full = MovieFine+repelem(MovieCoarse,CoarseSz/fineSz,CoarseSz/fineSz,1);


writerObj = VideoWriter('/Volumes/Lab/Users/bhaishahster/perturb.avi','Grayscale AVI');
writerObj.FrameRate=120;
open(writerObj);
writeVideo(writerObj,Movie_full);
close(writerObj);