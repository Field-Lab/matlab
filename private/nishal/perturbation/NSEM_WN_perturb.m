
%% Parameters
% Assume a grid of 320x640 


gridX=320;
gridY=640;

frameRate=120;
interval=1;
time=60; % in seconds
movieLen=time *frameRate/interval;

fineSz=2;

coarsemovieLen = time*frameRate/interval; 
CoarseSz = 2;

coarseframeRate = 1;

coarseScaleMax=0.5;

fineScale=0.24;
%% coarseMesh

load('/Volumes/Lab/Users/bhaishahster/NSbrownian_code_modifiable/matfiles/movie_chunk_1.mat');
movieChunkTime = size(movie,3)/(coarseframeRate*frameRate);
coarseRepeats = ceil(time / movieChunkTime);
CoarseMovieFine = repmat(movie,[1,1,coarseRepeats])*coarseScaleMax;
CoarseMovieFine = (permute (CoarseMovieFine,[2,1,3])/255) - 0.5;

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
    colorbar
    axis image
    
    subplot(1,2,2);
    imagesc(MovieFine(:,:,itime)');
    colormap gray
    axis image
    colorbar
    pause(1)
end