
%% Parameters
% Assume a grid of 320x640 


gridX=640;
gridY=320;

frameRate=120;
interval=1;
time=900; % in seconds
movieLen=time *frameRate/interval;

coarseinterval=1;
fineSz=4;

coarsemovieLen = time*frameRate/interval; 
CoarseSz = 8;

coarseType='checkerboard';
coarseframeRate = 1;

coarseScaleMax=0.36;
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
for itime=1:movieLen
     
            if(rem(itime,1000)==1)
                itime
            end
    for icntx=1:size(CoarseMovieFine,1)
        for icnty=1:size(CoarseMovieFine,2)
           
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

%Extend 
Movie_full = repmat(Movie_full,[2,2,1]);
CoarseMovieFine=repmat(CoarseMovieFine,[2,2,1]);

Movie_full=Movie_full(1:(gridX/fineSz) +16 ,1:(gridY/fineSz) +16 ,:);
CoarseMovieFine=CoarseMovieFine(1:(gridX/fineSz) +16 ,1:(gridY/fineSz) +16 ,:);
% 
% figure;
% for itime=1:movieLen
%     subplot(1,2,1);
%     imagesc(CoarseMovieFine(:,:,itime)');
%     colormap gray
%     axis image
%     
%     subplot(1,2,2);
%     imagesc(Movie_full(:,:,itime)');
%     colormap gray
%     axis image
%     pause(1)
% end

%% 
mov=uint8((Movie_full+0.5)*255);
save('/Volumes/Lab/Users/bhaishahster/NSbrownian_code_modifiable/Movie_full.mat','mov');
write_movie('/Volumes/Lab/Users/bhaishahster/NSbrownian_code_modifiable/Movie_full.mat','/Volumes/Lab/Users/bhaishahster/NSbrownian_code_modifiable/rawMovies/Movie_full.rawMovie',fineSz);


mov=uint8((CoarseMovieFine+0.5)*255);
save('/Volumes/Lab/Users/bhaishahster/NSbrownian_code_modifiable/CoarseMovieFine.mat','mov');
write_movie('/Volumes/Lab/Users/bhaishahster/NSbrownian_code_modifiable/CoarseMovieFine.mat','/Volumes/Lab/Users/bhaishahster/NSbrownian_code_modifiable/rawMovies/CoarseMovieFine.rawMovie',fineSz);