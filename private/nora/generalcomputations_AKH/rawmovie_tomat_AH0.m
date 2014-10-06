% DO NOT LOSE THIS !!  verions0  works but clunky 

% AK Heitman 2013-09-13
% Hopefully a slightly more robust/streamlined version of the code
% Reduces to the standard 80 by 40 8 pixel by 8 pixel version


clear; close all
%%
%rawmovie_dir = '/braid/snle/data/temporary_stimuli';
rawmovie_dir = '/braid/snle/scratch/Stimuli';
%rawmovie_dir = '/Akheitman/Stimuli';
decomp_dir = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/Decomposition';


% Lay down the parameters here
%{
rawmovie_name = 'eye-120-3_0-3600.rawMovie'; plots = true;
framesperblock = 120*30;
reduced_dim = [80 40];
blocks = 1:120;
params.moviename      = rawmovie_name;
params.reduced_dim    = reduced_dim; 
params.framesperblock = framesperblock;
%params.blocks         = blocks;
params.secondsperblock= 30;
params.refreshrate    = '120hz';
params.noteX          = 'unaltered movie.. in uint8 form.. goes from 0 to 256 I belive.. it is a linear scale';
params.noteXr         = 'downsampled movie  to a 80 by 40... equivalent to standard 8 pixel block used in white noise runs';
params.noteX_rescaled = 'everything shifted to a +- pt.5 scale for easier mathematical interpretation';
params.noteX          = 'all movies are index1 = frame, index2 = xval, index3 = yval';
FOI = 120:120:framesperblock;  plotsperblock =5; picsperplot = 6; subplot1 = 3; subplot2 = 2;
%}


%{
rawmovie_name = 'eye-long-v2.rawMovie'; plots = true;
framesperblock = 120*30;
reduced_dim = [80 40];
blocks = 1:120;
params.moviename      = rawmovie_name;
params.reduced_dim    = reduced_dim;
params.framesperblock = framesperblock;
%params.blocks         = blocks;
params.secondsperblock= 30;
params.refreshrate    = '120hz';
params.noteX          = 'unaltered movie.. in uint8 form.. goes from 0 to 256 I belive.. it is a linear scale';
params.noteXr         = 'downsampled movie  to a 80 by 40... equivalent to standard 8 pixel block used in white noise runs';
params.noteX_rescaled = 'everything shifted to a +- pt.5 scale for easier mathematical interpretation';
params.noteX          = 'all movies are index1 = frame, index2 = xval, index3 = yval';
FOI = 120:120:framesperblock;  plotsperblock =5; picsperplot = 6; subplot1 = 3; subplot2 = 2;

%}

rawmovie_name = 'eye-traces.rawMovie'; plots = true;
framesperblock = 4800;
reduced_dim = [80 40];
blocks = 1;
params.moviename      = rawmovie_name;
params.reduced_dim    = reduced_dim;
params.framesperblock = framesperblock;
%params.blocks         = blocks;
params.secondsperblock= 100;
params.refreshrate    = '120hz';
params.noteX          = 'unaltered movie.. in uint8 form.. goes from 0 to 256 I belive.. it is a linear scale';
params.noteXr         = 'downsampled movie  to a 80 by 40... equivalent to standard 8 pixel block used in white noise runs';
params.noteX_rescaled = 'everything shifted to a +- pt.5 scale for easier mathematical interpretation';
params.noteX          = 'all movies are index1 = frame, index2 = xval, index3 = yval';
FOI = 120:120:framesperblock;  plotsperblock = 40; picsperplot = 1; subplot1 = 1; subplot2 = 1;
%}

%{
rawmovie_name = 'eye-traces-2.rawMovie'; plots = true;
framesperblock = 12000;
reduced_dim = [80 40];
blocks = 1;
params.moviename      = rawmovie_name;
params.reduced_dim    = reduced_dim;
params.framesperblock = framesperblock;
%params.blocks         = blocks;
params.secondsperblock= 100;
params.refreshrate    = '120hz';
params.noteX          = 'unaltered movie.. in uint8 form.. goes from 0 to 256 I belive.. it is a linear scale';
params.noteXr         = 'downsampled movie  to a 80 by 40... equivalent to standard 8 pixel block used in white noise runs';
params.noteX_rescaled = 'everything shifted to a +- pt.5 scale for easier mathematical interpretation';
params.noteX          = 'all movies are index1 = frame, index2 = xval, index3 = yval';
FOI = 120:120:framesperblock;  plotsperblock = 25; picsperplot = 4; subplot1 = 2; subplot2 = 2;
%}

%{
rawmovie_name = 'eye-traces-2b.rawMovie'
plots = true
framesperblock = 12000;
reduced_dim = [80 40];
blocks = 1;
params.moviename      = rawmovie_name;
params.reduced_dim    = reduced_dim;
params.framesperblock = framesperblock;
params.blocks         = blocks;
params.secondsperblock= 30;
params.refreshrate    = '120hz';
params.noteX          = 'unaltered movie.. in uint8 form.. goes from 0 to 256 I belive.. it is a linear scale';
params.noteXr         = 'downsampled movie  to a 80 by 40... equivalent to standard 8 pixel block used in white noise runs';
params.noteX_rescaled = 'everything shifted to a +- pt.5 scale for easier mathematical interpretation';
params.noteX          = 'all movies are index1 = frame, index2 = xval, index3 = yval';
FOI = 120:120:framesperblock;  plotsperblock = 25; picsperplot = 4; subplot1 = 2; subplot2 = 2;
%}



%%
savedir = sprintf('%s/NSEM_directsaving_%s',decomp_dir,rawmovie_name)
if ~exist(savedir,'dir')
   unix(sprintf('mkdir %s',savedir));
end
eval(sprintf('save %s/Params.mat params', savedir) ) 

if plots
    plotdir = sprintf('%s/plots', savedir);
    if ~exist(plotdir, 'dir'), mkdir(plotdir); end
end


fid = fopen(sprintf('%s/%s',rawmovie_dir,rawmovie_name),'r'); % load movie
t = fscanf(fid,'%s',1);
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

%%

%%% NEEDS THIS TO RESTART THE READING
fid = fopen(sprintf('%s/%s',rawmovie_dir,rawmovie_name),'r'); 
%%% SKIP HEADER GO STRAIGHT TO GOOD STUFF %%%%%
fread(fid, header_size); % skip header
Xr = zeros(framesperblock,80,40); % used to be nan (not a number) instead of zeros

for k = blocks
   framenums = (k-1)*framesperblock+1:k*framesperblock;
   %  X = get_raw_movie2(f_mov,FR,[1,0,0]); % note: X is uint8
   
   %----------------
   
   if ~exist('X','var') % initialize
       X = zeros(framesperblock,width,height,'uint8');
   end
   
   
   % Roughly 30  seconds to load up in Bertha for a 30 second movie
   % Roughly 90 seconds on Alligator for a 30 second ovie
   % Loading up a Raw Movie 
   display(sprintf('Loading Full Movie of Block %d', k) );
   for i = 1:length(framenums)
         f = framenums(i);
         t = fread(fid,width*height*3,'ubit8');  % I think the 3 relates to RGB guns
         tt = reshape(t,3,width,height);
         X(i,:,:) = tt(1,:,:);
         %display(sprintf('pctdone of block = %d',i/max(framenums)))
   end
   
   %plotsperblock =5; picsperplot = 6;
   if plots
       for i_plot = 1:plotsperblock
           clf
           bigcount = (i_plot-1)*picsperplot;
           for i_pic = 1:picsperplot
               frame = FOI(bigcount + i_pic);
               subplot(subplot1,subplot2,i_pic)
               imagesc( (squeeze(X(frame,:,:))' )); colormap(gray);
           end
           orient tall
           eval(sprintf('print %s/Original_block_%d_image_%d.pdf -dpdf',plotdir,k,i_plot))
       end
   end
           
           
   
   %----------------
   
   % Takes roughly 60 seconds on Bertha    to load up a 30 second movie
   % Takes roughly 90 seconds on Alligator to load up a 30 second movie
   display(sprintf('Scaling Down to 80,40 Block %d', k) );
   for f = 1:framesperblock
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
   X_rescaled = ( (Xr-1) /255 ) - .5; 
  
   % Just 10 seconds to save 
   % 100 Megs per Block in this format
   eval(sprintf('save %s/eyemov%d X Xr X_rescaled',savedir,k))
end

fclose(fid);
