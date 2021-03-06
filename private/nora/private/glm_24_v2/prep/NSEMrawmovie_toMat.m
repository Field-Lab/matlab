% function      NSEMrawmovie_toMat(string_NSEMmoviename, frames_perblock, base_savedir, boolean_debug, boolean_NSEMFITplots)
%               2013-12-15  AKHeitman 
%
% usage:        Write .rawMovie into unit8 matrices 
%               Split into movieblocks of length frames_perblock
%               Plots all pictures to verify things are ok
%
% arguments:    string_NSEMmoviename
%
%
% calls:        movieinfo = NSEM_rawmovieinfo(moviename);
%
% outputs:      matfiles of the .rawMovie      
%
% paths:        run glmpath_18.m before everything
%
% oldversion    /generalcomputations_AKH/rawmovie_to_mat_AH 
%               will save here savedir =
%                 sprintf('%s/NSEM_%s/fullmovie_matfiles_%dframeblocks', base_savedir, moviename, framesperblock);
%%%%   NEEDS PROPER .RAWMOVIE CONVENTIONS   %%%%%%%%%%%  


% eg. calling: NSEMrawmovie_toMat('eye-210-3_0.rawMovie, framesperblock, base_savedir, boolean_NSEMFITplots)  

%% SETUP
function NSEMrawmovie_toMat(string_NSEMmoviename, framesperblock, base_savedir, boolean_NSEMFITplots)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  FOR DEBUGGING RUN FROM HERE %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
clear; close all
moviename = 'eye-long-v2';
framesperblock = 120*30;
boolean_debug = false;
boolean_NSEMFITplots = true;
base_savedir = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli';
rawmovie_otherdir = '/Akheitman/VanHateren/newrawmovie'
rawmovie_dir = '/braid/snle/scratch/Stimuli';
decomp_dir = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/Decomposition';
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%SETTING DIRECTORIES AND LOADING MOVIE INFO %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
moviename = string_NSEMmoviename;
savedir = sprintf('%s/NSEM_%s/fullmovie_matfiles_%dframeblocks', base_savedir, moviename, framesperblock);
plotdir = sprintf('%s/printedimages',savedir);
if ~exist(savedir,'dir'), mkdir(savedir); end
if boolean_NSEMFITplots && ~exist(plotdir, 'dir') 
    mkdir(plotdir); 
end
movieinfo              = NSEM_rawmovieinfo(moviename);
fullfile_location_name =  movieinfo.rawmoviefile ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% READ .RAWMOVIE PARAMS AND HEADER %%%%%%%%%%%%%
%%% NEEDS PROPER .RAWMOVIE CONVENTIONS %%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(fullfile_location_name,'r'); %get pointer for .rawMovie file
t = fscanf(fid,'%s',1);
% 
if ~isequal(t,'header-size')
   error('no header-size')
else
   header_size = str2double(fscanf(fid, '%s', 1));
end
height = [];
width = [];
while ( isempty(height) || isempty(width) )
   t = fscanf(fid,'%s',1);
   switch t
      case 'height'
         height = str2double(fscanf(fid,'%s',1));
      case 'width'
         width = str2double(fscanf(fid,'%s',1));
      otherwise
         fscanf(fid,'%s',1);
   end
end
if width ~= movieinfo.dimension(1)
    error('rawMovie dimesnion not matching information from /bookkeep/NSEM_rawmovieinfo')
end
if height ~= movieinfo.dimension(2)
    error('rawMovie dimesnion not matching information from /bookkeep/NSEM_rawmovieinfo')
end
% Determine number of blocks necessary 
movieblocks = floor(movieinfo.totalframes / framesperblock);
if movieblocks ~= (movieinfo.totalframes / framesperblock)
    display('frames per block is not perfect divisor.. Only first n blocks will be written.. last couple frames wil be dropped out of the matfiles')
end
fclose(fid);% CLOSE TO KEEP CLEAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% THIS IS NOT REDUNDANT ABSOLUTELY NECESSARY !!! %%%%%%%%%
%%% NEEDS THIS TO RESTART THE READING %%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(fullfile_location_name,'r'); 
fread(fid, header_size); % SKIP HEADER 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i_movieblock= 1:movieblocks
   
   framenums = (i_movieblock-1)*framesperblock+1:i_movieblock*framesperblock;   
   X = zeros(framesperblock,width,height,'uint8');

   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%% WRITING .RAWMOVIE TO UINT8-MATRIX %%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   % Roughly 60  seconds to load up in Bertha for a 30 second movie
   % Roughly 90 seconds on Alligator for a 30 second Movie
   display(sprintf('Loading Full Movie Block %d of %d', i_movieblock, movieblocks) );
   tic
   for i = 1:length(framenums)
         f = framenums(i);
         t = fread(fid,width*height*3,'ubit8');  % I think the 3 relates to RGB guns
         tt = reshape(t,3,width,height);
         X(i,:,:) = tt(1,:,:);
         %display(sprintf('pctdone of block = %d',i/max(framenums)))
   end
   toc
   %%% PACKAGE INTO MATFILE AND SAVE %%%%%%%%%%%%%%%
   movie.matrix = X;
   movie.movieblocknumber = i_movieblock;
   movie.framenumbers = framenums ;
   movie.note = 'Direct transcription from .rawMovie';
   movie.name = moviename;
   eval(sprintf('save %s/movieblock%d.mat movie', savedir, i_movieblock));
   %%% STANDARD PLOTTING FOR ALL NSEM MOVIES %%%%%%%%
   clf;
   subplot(2,1,1); imagesc( (squeeze(X(1,:,:)))' ); hold on; 
   colormap(gray); colorbar; xlabel('first frame'); hold off   
   subplot(2,1,2); imagesc( (squeeze(X(framesperblock,:,:)))' ); hold on;
   colormap(gray); colorbar; xlabel('last frame'); hold off   
   orient tall
   eval(sprintf('print %s/imageverify_block%d.pdf -dpdf', savedir, i_movieblock)) 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%% PLOTTING SPECFICALLY FOR NSEM-FIT MOVIES  %%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   if boolean_NSEMFITplots
       picsperplot = 6;
       plotsperblock = floor (framesperblock / (picsperplot* 120) );
       for i_plot = 1:plotsperblock
           clf
           bigcount = (i_plot-1)*picsperplot;
           for i_pic = 1:picsperplot
               frame = (i_plot-1)* (picsperplot* 120)  +(i_pic-1) * 120 +1;
               subplot(3,2,i_pic)
               imagesc( (squeeze(X(frame,:,:))' )); colormap(gray);
           end
           orient tall
           eval(sprintf('print %s/Original_block_%d_image_%d.pdf -dpdf',plotdir,i_movieblock,i_plot))
       end
       if plotsperblock ~= framesperblock / (picsperplot*120)
           clf
           remainder = floor( (framesperblock - (plotsperblock*picsperplot*120))/120);
           framestart = (plotsperblock*picsperplot*120) +1
           for i_pic = 1:remainder
               frame = framestart + 120*(i_pic-1);
               subplot(3,2,i_pic)
               imagesc( (squeeze(X(frame,:,:))' )); colormap(gray);
           end
           orient tall
           eval(sprintf('print %s/Original_block_%d_image_%d.pdf -dpdf',plotdir,i_movieblock,plotsperblock+1))
       end    
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
fclose(fid);
end
