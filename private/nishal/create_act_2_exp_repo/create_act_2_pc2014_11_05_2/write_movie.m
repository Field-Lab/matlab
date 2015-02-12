function write_movie(origin,destination,stixel)

% Use this to generate the rawMovie for use on the stim computers. Expects
% matfiles in a matfiles/movie_chunk_NUMBER.mat where NUMBER goes from 1 to
% number_images*segments_per_image
% this WILL overwrite a movie file if there is aleady one in the folder!

% movie settings

%stixel=10;
time_per_image=1;
movie_path = destination;%'/Volumes/Data/stimuli/movies/null_space/original_eyemove_stix10.rawMovie';

% size of image
w=1536;
h=1024;

% size of screen
screenw=320;% changed for Alex
screenh=320;

% number of pixels to be used for screen
wi=screenw/stixel;
hi=screenh/stixel;

% open the file to save the movie
[fid,msg] = fopen(movie_path, 'w');  

tic
 loaded_data=load(origin);%'/Volumes/Analysis/nishal/original/movie.mat');
 mov=loaded_data.mov;
 mov_len=size(mov,3);
 
% write header
hdr = sprintf('width\t%d\r\n', wi);
hdr = [hdr, sprintf('height\t%d\r\n', hi)];
hdr = [hdr, sprintf('frames-generated\t%d\r\n', mov_len)];
hdr = [hdr, sprintf('algorithm\teye_movement_matlab\r\n')];
hdr = [hdr, sprintf('images\tvanhateren\r\n')];
hdr = [hdr, sprintf('eye-movement\tbrownian\r\n')];
hdr = [hdr, sprintf('\r\n')];

fprintf(fid, sprintf('header-size\t%.10d\r\n', length(hdr)+24));
fprintf(fid, hdr);

%write frames

   
    for i=1:mov_len
        im=squeeze(mov(:,:,i));
        im = reshape(im,[1, numel(im)]);
        fwrite(fid, repmat(im, [3,1]), 'ubit8');
    end
    
fclose(fid);
toc
