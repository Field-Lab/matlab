% AKHeitman  2013-12-15
% Go from raw matrices of NSEM values to single mat file ready for GLM usage
% Organize into raster and fit blocks
% Hack to test whether we can gernalize to eye-longv2 ,  FEM900FF_longrast
%
%%%

%%%%%
%%%%  2013-12-16
% Absolutely needs to be redone to take scheming into account %
%
%  function( stimname, scheme)
%%%%
%%%%
%%%%


%Has problems outside of scheme A%%%

clear; close all;

NSEM_moviename = 'FEM900FF_longrast';
%NSEM_moviename = 'eye-long-v2';

movieinfo      = NSEM_rawmovieinfo(NSEM_moviename);
exp_nm = '2013-10-10-0'; fit_type = 'NSEM'; boolean_debug = false; map_type = 'mapPRJ';
[StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v18_split(exp_nm, boolean_debug, fit_type, map_type) 
SPars = StimulusPars.NSEM;

basedir = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_FEM900FF_longrast'
loaddir= '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli/NSEM_FEM900FF_longrast/fullmovie_matfiles_3600frameblocks';
fit_savename = sprintf('%s/blockedfitmovie_schemeB_8pixel.mat',basedir);
test_savename = sprintf('%s/testmovie_schemeB_8pix_Identity_8pix.mat',basedir);



%Scheme B
%
dim = [80,40,14400];
index{1} = [1:3600];
index{2} = [3601:7200];
index{3} = [7201:10800];
index{4} = [10801:14400];
%}



params.height = dim(2); 
params.width  = dim(1);
params.raster_frames = SPars.raster_frames;
params.fit_frames    = SPars.fit_frames;
params.rasterblocks        = SPars.StaticBlocks;
params.fitblocks           = SPars.NovelBlocks;
params.refresh_time        = SPars.Refresh;
params.refresh_time_note   = 'refreshtime is frame rate in msecs';


NSEMmovie.note0        = NSEM_moviename;   
NSEMmovie.note1        = 'Matrix is from 0 to 255 form of the movie, saved at uint8';
NSEMmovie.note2        = 'Move away from perturbative framework.. no light is 0, not a negative value';
NSEMmovie.note3        = 'This is downsampled to 8 by 8 CRT pixels (from 2 by 2 CRT) pixels';
NSEMmovie.note4        = 'fit and raster blocks are 120 seconds long';
NSEMmovie.params       = params;
NSEMmovie.fullParams   = SPars;
NSEMmovie.fitmovie.ind_to_block = SPars.NovelBlocks;







%%Scheme B  Fit Blocks
%{
for i_blk = SPars.NovelBlocks 
    display(sprintf('loading up %d', i_blk));
    matrix = zeros(dim);
    matfilenum = [i_blk , i_blk + 1];
    dim = [80,40,14400];
    index{1} = [1:3600];
    index{2} = [3601:7200];
    index{3} = [7201:10800];
    index{4} = [10801:14400];
    for i_half = 1:4
        
        clear Xr X movie       
        eval(sprintf('load %s/movieblock%d.mat movie',loaddir,matfilenum(i_half)));
        X  = movie.matrix;
        Xr = zeros(3600,80,40);
        
        for f = 1:3600
            tmp = double(squeeze(X(f,:,:)));
            for i = 1:80
                for j = 1:40
                    tmp2 = tmp(4*(i-1)+1:4*i,4*(j-1)+1:4*j);
                    Xr(f,i,j) = round(sum(tmp2(:))/16);  % make it integer valued
                end
            end
            if rem(f,360) == 0
                fprintf('now%d done\r',f)
            end
        end
       
        mov = Xr;
        mov2 = reshape(mov, [3600, 80*40]);
        mov3 = mov2';
        mov4 = reshape(mov3, [80 40 3600]);
        movf = (mov4);
        matrix(:,:,index{i_half}) = movf;
    end
    ind = i_blk/2;
    NSEMmovie.fitmovie.movie_byblock{ind}.matrix = uint8(matrix); 
end
eval(sprintf('save %s NSEMmovie',fit_savename));
%}

%%Scheme B Test Blk
for i_blk = 1 
    display(sprintf('loading up %d', i_blk));
    matrix = zeros(dim);
    matfilenum = [i_blk , i_blk + 1,i_blk +2, i_blk + 3];
    dim = [80,40,14400];
    index{1} = [1:3600];
    index{2} = [3601:7200];
    index{3} = [7201:10800];
    index{4} = [10801:14400];
    for i_half = 1:4
        
        clear Xr X movie       
        eval(sprintf('load %s/movieblock%d.mat movie',loaddir,matfilenum(i_half)));
        X  = movie.matrix;
        Xr = zeros(3600,80,40);
        
        for f = 1:3600
            tmp = double(squeeze(X(f,:,:)));
            for i = 1:80
                for j = 1:40
                    tmp2 = tmp(4*(i-1)+1:4*i,4*(j-1)+1:4*j);
                    Xr(f,i,j) = round(sum(tmp2(:))/16);  % make it integer valued
                end
            end
            if rem(f,360) == 0
                fprintf('now%d done\r',f)
            end
        end
       
        mov = Xr;
        mov2 = reshape(mov, [3600, 80*40]);
        mov3 = mov2';
        mov4 = reshape(mov3, [80 40 3600]);
        movf = (mov4);
        matrix(:,:,index{i_half}) = movf;
    end
    ind = i_blk/2;
    testmovie.matrix = uint8(matrix); 
    testmovie.params = params;
end
eval(sprintf('save %s testmovie',savename));
%}




%{
Scheme A
dim = [80,40,7200];
index{1} = [1:3600];
index{2} = [3601:7200];
for i_blk = SPars.NovelBlocks 
    display(sprintf('loading up %d', i_blk));
    matrix = zeros(dim);
    matfilenum = [i_blk , i_blk + 1];
    
    for i_half = 1:4
        
        clear Xr X movie       
        eval(sprintf('load %s/movieblock%d.mat movie',loaddir,matfilenum(i_half)));
        X  = movie.matrix;
        Xr = zeros(3600,80,40);
        
        for f = 1:3600
            tmp = double(squeeze(X(f,:,:)));
            for i = 1:80
                for j = 1:40
                    tmp2 = tmp(4*(i-1)+1:4*i,4*(j-1)+1:4*j);
                    Xr(f,i,j) = round(sum(tmp2(:))/16);  % make it integer valued
                end
            end
            if rem(f,360) == 0
                fprintf('now%d done\r',f)
            end
        end
       
        mov = Xr;
        mov2 = reshape(mov, [3600, 80*40]);
        mov3 = mov2';
        mov4 = reshape(mov3, [80 40 3600]);
        movf = (mov4);
        matrix(:,:,index{i_half}) = movf;
    end
    ind = i_blk/2;
    NSEMmovie.fitmovie.movie_byblock{ind}.matrix = uint8(matrix); 
end
eval(sprintf('save %s NSEMmovie',savename));
%}

%{
clear matrix
eval(sprintf('load %s/movieblock1.mat Xr',loaddir));
matrix = zeros(80,40,3600);
mov = Xr;
mov2 = reshape(mov, [3600, 80*40]);
mov3 = mov2';
mov4 = reshape(mov3, [80 40 3600]);
movf = uint8(mov4);
matrix = movf;
NSEMmovie.rastermovie.matrix = uint8(matrix);
%}



%{
Xr = zeros(dim);
 for f = 1:nfb
      tmp = double(squeeze(X(f,:,:)));
      for i = 1:80
         for j = 1:40
            tmp2 = tmp(4*(i-1)+1:4*i,4*(j-1)+1:4*j);
            Xr(f,i,j) = sum(tmp2(:))/16;
         end
      end
      if rem(f,360) == 0
         fprintf('now%d done\r',f)
      end
end
%}









%{
clear
exp_nm = '2012-08-09-3'; fit_type = 'BW'; boolean_debug = false; map_type = 'mapEI';
[StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v18_split(exp_nm, boolean_debug, fit_type, map_type) 
SPars = StimulusPars.NSEM;
%%

NSEM_moviename = 'eye-120-3_0-3600';
basedir = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli';
loaddir = sprintf('%s_oldstorage/NSEM_directsaving_%s',basedir, NSEM_moviename);
savedir = sprintf('%s/NSEM_%s',basedir,NSEM_moviename);

if ~exist(savedir, 'dir'), mkdir(savedir); end

dim = [80,40,7200];
index{1} = [1:3600];
index{2} = [3601:7200];

params.height = dim(2); 
params.width  = dim(1);
params.raster_frames = SPars.raster_frames;
params.fit_frames    = SPars.fit_frames;
params.rasterblocks        = SPars.StaticBlocks;
params.fitblocks           = SPars.NovelBlocks;
params.refresh_time     = SPars.Refresh;
params.refresh_time_note = 'refreshtime is frame rate in msecs';

NSEMmovie.note0        = NSEM_moviename;   
NSEMmovie.note1        = 'Matrix is from 0 to 255 form of the movie, saved at uint8';
NSEMmovie.note2        = 'Move away from perturbative framework.. no light is 0, not a negative value';
NSEMmovie.note3        = 'This is downsampled to 8 by 8 CRT pixels (from 2 by 2 CRT) pixels';
NSEMmovie.note4        = 'fit blocks are 60 seconds long, raster blocks are 30 seconds long';
NSEMmovie.params       = params;
NSEMmovie.fullParams   = SPars;
NSEMmovie.fitmovie.ind_to_block = SPars.NovelBlocks;
for i_blk = SPars.NovelBlocks 
    display(sprintf('loading up %d', i_blk));
    matrix = zeros(dim);
    matfilenum = [i_blk , i_blk + 1];
    
    for i_half = 1:2
        clear Xr mov mov2 mov3 mov4 movf
        eval(sprintf('load %s/eyemov%d.mat Xr',loaddir,matfilenum(i_half)));
        mov = Xr;
        mov2 = reshape(mov, [3600, 80*40]);
        mov3 = mov2';
        mov4 = reshape(mov3, [80 40 3600]);
        movf = (mov4);
        matrix(:,:,index{i_half}) = movf;
    end
    ind = i_blk/2;
    NSEMmovie.fitmovie.movie_byblock{ind}.matrix = uint8(matrix); 
end

% Raster movie
clear matrix
eval(sprintf('load %s/eyemov1.mat Xr',loaddir));
matrix = zeros(80,40,3600);
mov = Xr;
mov2 = reshape(mov, [3600, 80*40]);
mov3 = mov2';
mov4 = reshape(mov3, [80 40 3600]);
movf = uint8(mov4);
matrix = movf;
NSEMmovie.rastermovie.matrix = uint8(matrix);
eval(sprintf('save %s/NSEMmovie_8by8.mat NSEMmovie',savedir));
%}
