function unpack_rawmovie(moviepath, outputpath, framesperchunk)
%UNPACK_RAWMOVIE Unpacks raw movie files
%
%  UNPACK_RAWMOVIE(MOVIEPATH, OUTPUTPATH) unpacks a raw movie 
%  and stores all of its frames in individual .mat files in the specified 
%  output path. When called like this, the function assumes you are
%  unpacking a greyscale movie file.
%  By default, this function will store 120 frames per movie chunk. This
%  value can be optionally changed as follows.
%
%  UNPACK_RAWMOVIE(..., GREYSCALE) specifies if the movie is greyscale or 
%
%  UNPACK_RAWMOVIE(..., FRAMESPERCHUNK) will used the specified 
%  number of frames per mat file chunk.
%
%  Author: Georges Goetz - ggoetz@stanford.edu

if nargin == 2
    grayscale = true;
    framesperchunk = 120;
end
if nargin == 3
    framesperchunk = 120;
end

% Get the width, height and number of frames of the movie from the header
fid = fopen(moviepath,'r');
t = fscanf(fid,'%s',1);
if ~isequal(t,'header-size')
    error('no header-size')
else
    header_size = str2double(fscanf(fid, '%s', 1));
end

height = [];
width = [];
nframes = [];
while ( isempty(height) || isempty(width) || isempty(nframes) )
    t = fscanf(fid,'%s',1);
    switch t
        case 'height'
            height = str2double(fscanf(fid, '%s', 1));
        case 'width'
            width = str2double(fscanf(fid, '%s', 1));
        case 'frames-generated'
            nframes = str2double(fscanf(fid, '%s', 1));
        otherwise
            fscanf(fid,'%s',1);
    end
end

% Figure out the file name 
if exist(outputpath) %#ok<EXIST>
    assert(isdir(outputpath))
else
    mkdir(outputpath)
end
basefilename = split(outputpath, filesep);
basefilename = basefilename{end};
filename = sprintf('%s_%c0%dd.mat', ...
    basefilename, hex2dec('25'), ...
    length(num2str(ceil(nframes/framesperchunk))));

% Now we can read the movie 
fid = fopen(moviepath, 'r');
fread(fid, header_size); % skip header

% Create a crude progress bar
fprintf('Converting movie to mat files\n');
fprintf([repmat('.',1,80) '\n']);
ndots = 0;

% Loading up the Raw Movie
frames = cell(framesperchunk, 1);
for i = 1:nframes
    % Progress bar
    if (i/nframes)*80 > ndots
        fprintf('.');
        ndots = ndots + 1;
    end
    
     % Load
    t = fread(fid,width*height*3,'ubit8');  % 3 RGB guns
    tt = reshape(t,3,width,height);
    tt = permute(tt, [2, 3, 1]);
    
    if mod(i, framesperchunk) == 0
        if greyscale
            frames{framesperchunk} = squeeze(tt(:,:,1));
        else
            frames{framesperchunk} = tt;
        end
        save(fullfile(outputpath, sprintf(filename, i/framesperchunk)), 'frames');
        frames = cell(framesperchunk, 1);
    else
        if greyscale
            frames{framesperchunk} = squeeze(tt(:,:,1));
        else
            frames{framesperchunk} = tt;
        end
    end
end
% There could be a few last frames that haven't been saved yet
if mod(i, framesperchunk) ~= 0
    frames(cellfun(@isempty,frames)) = [];
    save(fullfile(outputpath, sprintf(filename, ceil(i/framesperchunk))), 'frames');
end
% Done!
fprintf('\nConversion done.\n');

fclose(fid);

end % unpack_rawmovie