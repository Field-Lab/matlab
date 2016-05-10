% Nora Brackbill 08-12-2014
% Use this to generate the rawMovie for use on the stim computers. Expects
% matfiles in a matfiles/movie_chunk_NUMBER.mat where NUMBER goes from 1 to
% number_images*segments_per_image
% this WILL overwrite a movie file if there is aleady one there!

function generate_rawmovie_from_mat(folder)
% load movie param
% where to save raw movie
movie_path = [folder '/newrawmovie/'];

% size of screen
screenw=640;
screenh=320;

% number of pixels to be used for screen
wi=screenw/2;
hi=screenh/2;

% open the file to save the movie
[fid,msg] = fopen(movie_path, 'w');  

tic

% write header
hdr = sprintf('width\t%d\r\n', wi);
hdr = [hdr, sprintf('height\t%d\r\n', hi)];
hdr = [hdr, sprintf('frames-generated\t%d\r\n', params.rawmovie_number_images*120)];
hdr = [hdr, sprintf('algorithm\teye_movement_matlab\r\n')];
hdr = [hdr, sprintf('images\tvanhateren\r\n')];
hdr = [hdr, sprintf('eye-movement\tbrownian\r\n')];
hdr = [hdr, sprintf('\r\n')];

fprintf(fid, sprintf('header-size\t%.10d\r\n', length(hdr)+24));
fprintf(fid, hdr);

%write frames
all_files_read = 0;
while ~all_files_read
    try
        load([folder '/matfiles/movie_chunk_' num2str(j)]);
        display(sprintf('imagenumber %d', j))
        for i=1:size(movie_chunk,3)
            im=squeeze(movie_chunk(:,:,i));
            im = reshape(im,[1, numel(im)]);
            fwrite(fid, repmat(im, [3,1]), 'ubit8');
        end
        j = j+1;
    catch
        all_files_read = 1; 
    end
end
fclose(fid);
toc
end
