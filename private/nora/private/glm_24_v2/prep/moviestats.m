% Hack code AkHeitman  2014-02-04 final
% Produces correlation strucure, distribution histogram /mu etc.
% Can be used to derive optimal linear filters.. avg image  etc.
 
%%
% Movie Analysis
clear; close all

dBug = false;
%stimname= 'eye-120-3_0-3600'; cmodel = '8pix_Identity_8pix'; blocks = 59; fpb = 7200;  scheme = 'schemeA';
%stimname= 'eye-120-3_0-3600'; cmodel = '8pix_Model1_1e4_8pix'; blocks = 59; fpb = 7200;  scheme = 'schemeA';
%stimname= 'eye-120-3_0-3600'; cmodel = '8pix_Model1_1e6_8pix'; blocks = 59; fpb = 7200;  scheme = 'schemeA';
stimname= 'eye-120-3_0-3600'; cmodel = '8pix_Model1_1e5_8pix'; blocks = 59; fpb = 7200;  scheme = 'schemeA';
%stimname= 'eye-long-v2'; cmodel = '8pix_Identity_8pix'; blocks = 59; fpb = 7200;  scheme = 'schemeA';
%stimname= 'eye-long-v2'; cmodel = '8pix_Model1_1e4_8pix'; blocks = 59; fpb = 7200;  scheme = 'schemeA';
%stimname = 'FEM900FF_longrast'; cmodel = '8pix_Model1_1e4_8pix'; blocks = 27; fpb = 14400; scheme = 'schemeB';
%stimname = 'FEM900FF_longrast'; cmodel = '8pix_Identity_8pix'; blocks = 27; fpb = 14400; scheme = 'schemeB';
%fpb = 14400; scheme = 'schemeB';
stimdir = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli';
moviedir = sprintf('%s/NSEM_%s', stimdir,stimname);
moviefile = sprintf('%s/fitmovie_%s_%s.mat', moviedir, scheme,cmodel); 
eval(sprintf('load %s', moviefile))
dim1 = 80; dim2 = 40;


display(sprintf('~~~~~Working on %s ~~~~', moviefile))
if dBug, blocks = 3; end
movie = (zeros(dim1,dim2,fpb * blocks));
for i_blk = 1 : blocks
    framea = (i_blk-1)*fpb + 1;
    frameb = i_blk * fpb;
    display(sprintf('Loading up frame %d, %d pct done', frameb, round(100*(frameb/(fpb*blocks) ))) );
	movie_singleform0 = double(NSEMmovie.fitmovie.movie_byblock{i_blk}.matrix);
    movie_singleform  = squeeze(movie_singleform0(:,:,:));
    movie_singleform = (1/255) * movie_singleform;
    movie(:,:,framea:frameb) = movie_singleform;
end
clear movie_singleform movie_singleform0 framea frameb i_blk



%%
%{
singlemovie = NSEMmovie.fitmovie.movie_byblock{1}.matrix;

movstruct(360).colormap = []; 

for f_num = 1:360
    imagemat = squeeze(singlemovie(:,:,f_num));
    cdatamat = cat(3,imagemat,imagemat);
    cdatamat = cat(3,cdatamat,imagemat);
    
    F(f_num).cdata = cdatamat;
    F(f_num).colormap = [];
end


figure;
tic
for f_num = 1:7200
    imagemat = squeeze(singlemovie(:,:,f_num));
    imagesc(imagemat'); colormap gray; 
    
    movstruct(f_num) = getframe;
end
toc
figure;
movie(movstruct,2);
%}

%%  Raw stats for the 0-1 scaled movie
inputstats.stimulus       = stimname;
inputstats.cmodel         = cmodel;
inputstats.dim            = [dim1, dim2];

movie_spt = reshape(movie, [dim1*dim2, fpb * blocks]);
frames = 120:120:(fpb*blocks);
npix = length(frames); 
pix_spt = zeros([dim1*dim2], npix);
for i_pix = 1:npix
    i_frame = frames(i_pix);
    pix_spt(:,i_pix) = movie_spt(:,i_frame);
end
inputstats.lastframe_persec  = pix_spt;




meanintensities_perpix = mean(movie_spt,2);
inputstats.full_avgIperpix = meanintensities_perpix;
inputstats.mu_avgIperpix = mean(meanintensities_perpix);
inputstats.std_avgIperpix = std(meanintensities_perpix);

display('Computing histogram , takes ')  % ~1minute bertha
tic
if strcmp(cmodel, '8pix_Identity_8pix')
    hist_8bit   =  hist(pix_spt(:) , 256) ; 
    inputstats.hist_8bit   =  (1 / (npix * dim1*dim2) ) * hist_8bit; 
else
    hist_8bit   =  hist(movie_spt(:) , 256) ; 
    inputstats.hist_8bit   =  (1 / (fpb*blocks * dim1*dim2) ) * hist_8bit;
end
toc

%% Correlation stats .. subtract out mean
meanzeromov  = movie_spt - repmat(meanintensities_perpix,1,fpb * blocks);
meanzeropics = pix_spt - repmat( (mean(pix_spt,2)) ,1,npix);

pixel_subset0 = round(dim1*rand(1,dim2));
pixel_subset  = pixel_subset0 + [0:dim1: (dim2-1)*dim1];
netsignal = zeros(1,120);
netsignal2 = zeros(1,120);
for i_pix = 1 : length(pixel_subset);
    x_sig = meanzeromov(pixel_subset(i_pix), :);
    cc = xcorr(x_sig ,150, 'biased');
    cc2 = xcorr(x_sig ,150, 'coeff');
    decay  = cc(151:270);
    decay2  = cc2(151:270);
   % temp_autocorr{i_pix}.autocorr  = decay;
    %temp_autocorr{i_pix}.autocorr2  = decay2;
    netsignal = netsignal + decay;
    netsignal2 = netsignal2 + decay2;
end
avgsignal = netsignal / length(pixel_subset);
avgsignal2 = netsignal2 / length(pixel_subset);
inputstats.temporalcorrelation_coeff = avgsignal2;  
inputstats.temporalcovariance        = avgsignal; 


display('Working on Spatial Covariance');
tic
ccc = cov(meanzeromov');
ddd = corrcoef(meanzeromov');
toc

inputstats.spatialcovariance = ccc;
inputstats.spatialcorrelation = ddd;
inputstats.inv_spatialcovariance = inv(inputstats.spatialcovariance); 


figure;
clims = [0 1];
subplot(3,2,1); plot( linspace(0,1,256) ,  inputstats.hist_8bit ); hold on
plot(inputstats.mu_avgIperpix, 0 , 'r*'); title('Normed Distribution'); ylabel('8bit bin');
subplot(3,2,2); plot( linspace(0,1,120) , inputstats.temporalcorrelation_coeff); hold on; xlabel('seconds'); title('Avg Pixel Auto-Corr');
subplot(3,2,3); imagesc(inputstats.spatialcorrelation , clims); colorbar; title('spatial corrcoef');
subplot(3,2,4);  imagesc( reshape(inputstats.spatialcorrelation(801,:),[dim1,dim2]), clims );  colorbar;  title('pt801');
subplot(3,2,5);  imagesc( reshape(inputstats.spatialcorrelation(1620,:),[dim1,dim2]),clims  ); colorbar;  title('pt1620');
subplot(3,2,6);  imagesc( reshape(inputstats.spatialcorrelation(3200,:),[dim1,dim2]),clims  ); colorbar;  title('pt3200');
orient tall
filename = sprintf('%s/fitmovie_structure_%s', moviedir,  cmodel);

if ~dBug
eval(sprintf('print -dpdf %s.pdf', filename))
eval(sprintf('save %s/inputstats_%s.mat inputstats', moviedir, cmodel) );
end


if dBug
eval(sprintf('print -dpdf %s_DEBUG.pdf', filename))
eval(sprintf('save %s/inputstats_%s_DEBUG.mat inputstats', moviedir, cmodel) );
end



%%


%%

