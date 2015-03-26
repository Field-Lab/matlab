%% Movie info
% AKHeitman 2014-01-26  hacked code
% RUN MUCH FASTER ON BERTHA!!!!



clear; close all;

%origmovie = 'eye-120-3_0-3600'; origscheme = 'schemeA';
origmovie = 'FEM900FF_longrast'; origscheme = 'schemeB';
%origmovie = 'eye-long-v2'; origscheme = 'schemeA';
moviedir = sprintf('/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_%s', origmovie);

eval(sprintf('load %s/testmovie_%s_8pix_Identity_8pix.mat', moviedir,origscheme)) ;
eval(sprintf('load %s/fitmovie_%s_8pix_Identity_8pix.mat', moviedir,origscheme)) ;

%maxRstar_sec = 1000000; model = 'Model1'; model_mode ='8pix_Model1_1e6_8pix'
%maxRstar_sec = 100000; model = 'Model1'; model_mode ='8pix_Model1_1e5_8pix'
maxRstar_sec = 10000;  model = 'Model1'; model_mode ='8pix_Model1_1e4_8pix';
savedir = sprintf('%s/%s_%s', moviedir,origscheme, model_mode) ;
if ~exist(savedir, 'dir'), mkdir(savedir); end

%% Run Model on each block
%{
%blocks = [1:30];
%blocks = [16:30];
%blocks = [31:45];
%blocks = [46:59];
%testblock = true; blocks = [1:20];
%testbock = false; blocks = [21:40]
%testbock = false;% blocks = [41:59]
testblock = true; blocks = [];
%%% CMP Cone Model Params 
binsperframe = 8  % roughly msec bins 
framedur = (1/120);
max_rstar_sec_mult = maxRstar_sec/255;
% tries to simulate having a half sec of gray screen (127 on the 255 scale) for the time preceeding each block
padframes = 60;
padval = 127;

CMP.padval = padval; CMP.padframes = padframes; CMP.max_rstar_sec_mult = max_rstar_sec_mult; 
CMP.binsperframe = binsperframe; CMP.framedur = framedur;
% Look at the 


%%% Doing the test movie %%%
if exist('testblock') && testblock
    display('--------- Working on the testblock -----------')
    moviematrix0 = testmovie.matrix;
    frontpad = padval*ones(size(moviematrix0,1),size(moviematrix0,2), padframes);
    prepad = double(moviematrix0);
    moviematrix = cat(3,frontpad, prepad);
    movie_rstar_frame = framedur * max_rstar_sec_mult * moviematrix;
    checkfunction.check = true;
    checkfunction.plotdir = sprintf('%s/checkplots', savedir);
    if ~exist(checkfunction.plotdir, 'dir'), mkdir(checkfunction.plotdir); end

    checkfunction.plotblockname = sprintf('rasterblock');

    model_name  = 'model1';
     % clear moviematrix padframes padval frontpad prepad max_rstar_sec_mult moviematrix0 
    pAmpmovie = runconemodel(model_name, movie_rstar_frame, framedur, binsperframe, checkfunction);
    cmodeloutput.pAmp         = pAmpmovie(:,:,padframes+1:end);
    cmodeloutput.block = 'rasterblock';
    cmodeloutput.movieblock   = 1;
    cmodeloutput.cmodel       = model_mode;
    cmodeloutput.params = CMP;
    eval(sprintf('save %s/cmodel_rasterblock.mat cmodeloutput', savedir));
end



%
for i_blk = blocks
   % display(sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'))
    display(sprintf('------------- Working on Block %d -------------', i_blk))   
   % display(sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'))
    moviematrix0  = NSEMmovie.fitmovie.movie_byblock{i_blk}.matrix;
    movieblocknum = NSEMmovie.fitmovie.ind_to_block(i_blk); 
    

    
    
    frontpad = padval*ones(size(moviematrix0,1),size(moviematrix0,2), padframes);
    prepad = double(moviematrix0);
    moviematrix = cat(3,frontpad, prepad);
    movie_rstar_frame = framedur * max_rstar_sec_mult * moviematrix;
    checkfunction.check = true;
    checkfunction.plotdir = sprintf('%s/checkplots', savedir);
    if ~exist(checkfunction.plotdir, 'dir'), mkdir(checkfunction.plotdir); end
    checkfunction.plotblockname = sprintf('fitblock%d', i_blk);
    model_name  = 'model1';
    pAmpmovie = runconemodel(model_name, movie_rstar_frame, framedur, binsperframe, checkfunction);
    
    
    cmodeloutput.pAmp         = pAmpmovie(:,:,padframes+1:end);
    cmodeloutput.fitblock     = i_blk;
    cmodeloutput.movieblock   = movieblocknum;
    cmodeloutput.cmodel       = model_mode;
    cmodeloutput.params = CMP;
    
    
    eval(sprintf('save %s/cmodel_fitblock_%d.mat cmodeloutput', savedir, i_blk));
end
%}






%%%
%% Collect block. .. find normalization constants
%
%
%
%
loaddir = savedir;
display(sprintf('Working on %s with model %s', origmovie, model_mode) );
NSEMmovie.note0 = model_mode;
NSEMmovie.note1 = 'original cmodel output is in units of pAmp';
NSEMmovie.note2a = 'first normalize to 0 1 to make unitless';
NSEMmovie.note2b = 'rescaled to 0 255 to save as uint8';
NSEMmovie.note2c = 'rescaling purely to save space ..';
NSEMmovie.note3 = 'orginal movie eye-120-3_0-3600.rawMovie';
NSEMmovie.params.width = 80;
NSEMmovie.params.height = 40;

if strcmp(origscheme, 'schemeA')
    ind_to_block = 2:2:118;
    NSEMmovie.params.fit_frames = 7200;
end

if strcmp(origscheme, 'schemeB')
    ind_to_block = 2:2:54;
    NSEMmovie.params.fit_frames = 14400;
end

NSEMmovie.normalization.pAmpmin  = zeros(size(ind_to_block));
NSEMmovie.normalization.pAmpmax  = zeros(size(ind_to_block));
NSEMmovie.normalization.pAmpmean = zeros(size(ind_to_block));
NSEMmovie.normalization.pAmpstd = zeros(size(ind_to_block));



% load test block too
display(sprintf('loading raster movie to find normalization constants'));
eval(sprintf('load %s/cmodel_rasterblock', loaddir));
mat0  = cmodeloutput.pAmp;

NSEMmovie.normalization.rast_pAmpmin   = min(mat0(:));
NSEMmovie.normalization.rast_pAmpmax   = max(mat0(:));
NSEMmovie.normalization.rast_pAmpmean  = mean(mat0(:));
NSEMmovie.normalization.rast_pAmpstd   = std(mat0(:));

% 1 Load up all blocks to find the   
for i_blk = 1:length(ind_to_block)
    display(sprintf('loading fit movie number %d out of 59 for normalization constants', i_blk));
    eval(sprintf('load %s/cmodel_fitblock_%d', loaddir,i_blk));
    mat0  = cmodeloutput.pAmp;
    NSEMmovie.normalization.pAmpmin(i_blk)  = min(mat0(:));
    NSEMmovie.normalization.pAmpmax(i_blk)  = max(mat0(:));
    NSEMmovie.normalization.pAmpmean(i_blk) = mean(mat0(:));
    NSEMmovie.normalization.pAmpstd(i_blk) = std(mat0(:));
   % clear cmodeloutput
end

absmin = min( min(NSEMmovie.normalization.pAmpmin) , NSEMmovie.normalization.rast_pAmpmin )
absmax = max( max(NSEMmovie.normalization.pAmpmax) , NSEMmovie.normalization.rast_pAmpmax )

NSEMmovie.normalization.minval    = absmin;
NSEMmovie.normalization.maxval    = absmax;
NSEMmovie.normalization.meanval   = mean(NSEMmovie.normalization.pAmpmean);


display(sprintf('loading raster movie to for unit rescaling'));
eval(sprintf('load %s/cmodel_rasterblock', loaddir));
mat0  =  cmodeloutput.pAmp;
mat1  =  mat0 - absmin;
mat2  = ((mat1/(absmax-absmin) )  ) ;
testmovie.normalization = NSEMmovie.normalization;
testmovie.params = NSEMmovie.params;
testmovie.matrix = mat2;  
eval(sprintf('save %s/testmovie_%s_%s.mat testmovie', moviedir, origscheme, model_mode));


for i_blk = 1:length(ind_to_block)
    display(sprintf('loading fit movie number %d out of 59 for rescaling and uint8 storage', i_blk));
    eval(sprintf('load %s/cmodel_fitblock_%d', loaddir,i_blk));
    mat0  = cmodeloutput.pAmp;
    mat1  = mat0 - absmin;
    mat2  = round(    255* (mat1/(absmax-absmin) )  ) ;
    NSEMmovie.fitmovie.movie_byblock{i_blk}.matrix = uint8(mat2);
end
eval(sprintf('save %s/fitmovie_%s_%s.mat NSEMmovie', moviedir,origscheme, model_mode));



%}


%clear; close all   
%{
moviedir = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_eye-120-3_0-3600';
testmovie.note0 = '8pix_model1_8pix';
testmovie.note1 = 'original cmodel output is in units of pAmp';
testmovie.note2 = 'rescaled to 0 255 to save as uint8';
testmovie.note3 = 'orginal movie eye-120-3_0-3600.rawMovie';


eval(sprintf('load %s/schemeA_8pix_model1_8pix/cmodel_rasterblock', moviedir));
mat0  = cmodeloutput.pAmp;
mat1  = mat0 - min(mat0(:));
mat2  = round(255* (mat1 / max(mat1(:))) );   
testmovie.matrix = uint8(mat2);
eval(sprintf('save %s/testmovie_schemeA_8pix_Model1_8pix.mat testmovie', moviedir));   
%}





