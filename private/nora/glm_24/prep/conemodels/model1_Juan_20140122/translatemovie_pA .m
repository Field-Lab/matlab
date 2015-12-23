% ahackin around 
load blockedmovie_schemeA_8pixel.mat

rawmatrix = double(NSEMmovie.fitmovie.movie_byblock{1}.matrix);
frames = size(rawmatrix,3);
Rstar_sec = 10000 / 255;
framedur = (1/120);
time0 = framedur : framedur : (frames * framedur);
xval =1 ; yval= 1;
Rstar_frame = framedur * Rstar_sec * rawmatrix;
pixmovie0 = squeeze(Rstar_frame(xval,yval,:))';


%%
binsperframe = 4;
bpf = binsperframe;
time = framedur/bpf : framedur/bpf : frames * framedur;




pixmovie = (1/binsperframe)*reshape( repmat(pixmovie0, bpf,1) , bpf*length(pixmovie0),1) ;
pixmovie = pixmovie';figure; 
figure; subplot(2,1,1); plot(time0,pixmovie0); subplot(2,1,2); plot(time, pixmovie); hold on; title('Check upsampling worked')   



%%
display(sprintf('time to run a 60 second pixel at %d bins per frame',binsperframe)) 
tic
afit=-riekemodelstim_slowadapt([DimFlash.coeffs(5) EyeMovements.slowcacoeff],time,pixmovie, DimFlash.coeffs(1:7));
toc
%%

stimnormed = (pixmovie - min(pixmovie)) / max(pixmovie - min(pixmovie));
fitnormed  =  (afit - min(afit)) / max(afit - min(afit));


origstim_normed = (pixmovie0 - min(pixmovie0)) / max(pixmovie0 - min(pixmovie0));
fitnormed_ds = interp1(time,fitnormed,time0); 


figure;
subplot(4,1,1); plot( time ,  pixmovie, 'b'); ylabel('R*perbin')
subplot(4,1,2); plot( time , afit, 'r'); ylabel('pA')
subplot(4,1,3); hold on; plot(time ,stimnormed,'b'); plot(time,fitnormed, 'r')
subplot(4,1,4); hold on; plot(time0 ,origstim_normed,'b'); plot(time0,fitnormed_ds, 'r')