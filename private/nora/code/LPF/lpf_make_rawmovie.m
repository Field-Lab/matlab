% Nora Brackbill 08-12-2014
% Use this to generate the rawMovie for use on the stim computers. Expects
% matfiles in a matfiles/movie_chunk_NUMBER.mat where NUMBER goes from 1 to
% seconds
% this WILL overwrite a movie file if there is aleady one there!

function lpf_make_rawmovie(matfiles_folder, movie_path, seconds, dim)

% number of pixels to be used for screen
wi=dim(1);
hi=dim(2);

% open the file to save the movie
[fid,msg] = fopen(movie_path, 'w');  

tic

% write header
hdr = sprintf('width\t%d\r\n', wi);
hdr = [hdr, sprintf('height\t%d\r\n', hi)];
hdr = [hdr, sprintf('frames-generated\t%d\r\n', seconds*120)];
hdr = [hdr, sprintf('algorithm\tLPF_v2\r\n')];
% hdr = [hdr, sprintf('images\tvanhateren\r\n')];
% hdr = [hdr, sprintf('eye-movement\tbrownian\r\n')];
hdr = [hdr, sprintf('\r\n')];

fprintf(fid, sprintf('header-size\t%.10d\r\n', length(hdr)+24));
fprintf(fid, hdr);

%write frames
for j=1:seconds
    load([matfiles_folder 'movie_chunk_' num2str(j)]);
    display(sprintf('second %d', j))
    for i=1:120
        im=squeeze(movie(:,:,i));
        im = reshape(im,[1, numel(im)]);
        fwrite(fid, repmat(im, [3,1]), 'ubit8');
    end
end
fclose(fid);
toc
end