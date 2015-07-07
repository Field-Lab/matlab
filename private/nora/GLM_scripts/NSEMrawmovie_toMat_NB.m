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
function NSEMrawmovie_toMat_NB(moviename, raster_seconds, fit_seconds, movieblocks)


basedir = '/Volumes/Data/Stimuli/movies/eye-movement/old_movies/';
savedir = ['/Volumes/Lab/Users/Nora/NSEM_Movies/' moviename];
if ~exist(savedir,'dir'), mkdir(savedir); end
fullfile_location_name =  [basedir moviename '.rawMovie'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FOR TESTING
% fullfile_location_name =  '/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian_6000/NSbrownian_6000_A_025.rawMovie';
raster_frames = raster_seconds * 120;
fit_frames = fit_seconds * 120;

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
fclose(fid);% CLOSE TO KEEP CLEAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% THIS IS NOT REDUNDANT ABSOLUTELY NECESSARY !!! %%%%%%%%%
%%% NEEDS THIS TO RESTART THE READING %%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(fullfile_location_name,'r'); 
fread(fid, header_size); % SKIP HEADER 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First get the raster

framenums = 1:raster_frames;
X = zeros(raster_frames,width,height,'uint8');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% WRITING .RAWMOVIE TO UINT8-MATRIX %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Roughly 60  seconds to load up in Bertha for a 30 second movie
% Roughly 90 seconds on Alligator for a 30 second Movie
display(sprintf('Loading Raster Block'));
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
movie.movieblocknumber = 0;
movie.framenumbers = framenums ;
movie.note = 'Direct transcription from .rawMovie';
movie.name = moviename;
eval(sprintf('save %s/testmovie.mat movie', savedir));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i_movieblock= 1:movieblocks
   framenums = raster_frames+((i_movieblock-1)*fit_frames+1:i_movieblock*fit_frames);   
   X = zeros(fit_frames,width,height,'uint8');

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

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
fclose(fid);
end
