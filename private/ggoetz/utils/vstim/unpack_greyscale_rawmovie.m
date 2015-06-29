function unpack_greyscale_rawmovie(moviepath, outputpath, framesperchunk)
%UNPACK_GREYSCALE_RAWMOVIE Unpacks greyscale raw movie files
%
%  UNPACK_GREYSCALE_RAWMOVIE(MOVIEPATH, OUTPUTPATH) unpacks a raw movie 
%  and stores all of its frames in individual .mat files in the specified 
%  output path. 
%  By default, this function will store 120 frames per movie chunk. This
%  value can be optionally changed as follows.
%
%  UNPACK_GREYSCALE_RAWMOVIE(..., FRAMESPERCHUNK) will used the specified 
%  number of frames per mat file chunk.
%
%  WARNING: this function should not be used for color raw movies. It will
%  return incorrect results if used for this purpose.

if nargin == 2
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
    frames{i} = squeeze(tt(1,:,:));
    
    if mod(i, framesperchunk) == 0
        save(fullfile(outputpath, sprintf(filename, i/framesperchunk)), 'frames');
        frames = cell(framesperchunk, 1);
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

end % unpack_greyscale_rawmovie