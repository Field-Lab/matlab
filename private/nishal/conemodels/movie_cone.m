function pAmpmovie = movie_cone(qq,interval)

moviematrix0 = (qq+0.5)*255;

binsperframe = 8;  %roughly msec bins
%interval=2;
framedur = (interval*1/120);
max_rstar_sec_mult = 1000/255;


% Front Pad by 1 s
padframes = 60;
padval = 127;
frontpad = padval*ones(size(moviematrix0,1),size(moviematrix0,2), padframes);
prepad = double(moviematrix0);
moviematrix = cat(3,frontpad, prepad);
movie_rstar_frame = framedur * max_rstar_sec_mult * moviematrix;
checkfunction.check = true;
checkfunction.plotdir = '/Volumes/Analysis/nishal/conemodel/';
checkfunction.plotblockname = sprintf('fitblock%d', 1);
model_name  = 'model1';
clear moviematrix padframes padval frontpad prepad max_rstar_sec_mult moviematrix0 


pAmpmovie = runconemodel(model_name, movie_rstar_frame, framedur, binsperframe, checkfunction);

figure;
for itime=1:size(pAmpmovie,3)
    itime
imagesc(pAmpmovie(:,:,itime));
colormap gray
pause(1/60)
end


pixX=30;
pixY=40;
timeindices = 1:size(pAmpmovie,3)-60;
conesig = squeeze(pAmpmovie(pixY,pixY,61:end));

conesignormed = (conesig - mean(conesig))/max(abs(conesig));
plotyy(timeindices,squeeze(qq(pixX,pixY,timeindices)),timeindices,conesig);
hold on;
plot(squeeze(qq(pixX,pixY,timeindices)))

pAmpmovie=pAmpmovie(:,:,61:end);
