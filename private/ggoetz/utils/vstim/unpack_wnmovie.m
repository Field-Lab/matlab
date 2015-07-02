function unpack_wnmovie(moviepath, globalspath, outputpath, nframes, framesperchunk)
%UNPACK_WNMOVIE Unpacks a white noise movie
%
%  UNPACK_WNMOVIE(MOVIEPATH, GLOBALSPATH, NFRAMES, OUTPUTPATH) unpacks a 
%  raw movie and stores all of its frames in individual .mat files in the 
%  specified output path. 
%  By default, this function will store 120 frames per movie chunk. This
%  value can be optionally changed as follows.
%  MOVIEPATH should be the path to a Vision .movie file and GLOBALSPATH
%  the path to a vision globals file. This function generates NFRAMES 
%  frames of the movie
%
%  UNPACK_GREYSCALE_RAWMOVIE(..., FRAMESPERCHUNK) will used the specified 
%  number of frames per mat file chunk.
%

if nargin == 4
    framesperchunk = 120;
end

% Create an instance of a vision movie file
wnmovie = edu.ucsc.neurobiology.vision.stimulus.WhiteNoiseMovie(moviepath, globalspath);

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
    
    % Load. Note: frames are 0-indexed (java convention)
    framedata = double(reshape(wnmovie.getFrame(i-1).getBuffer(), ...
        wnmovie.getHeight(), wnmovie.getWidth(), 3));
    
    if mod(i, framesperchunk) == 0
        frames{framesperchunk} = framedata;
        save(fullfile(outputpath, sprintf(filename, i/framesperchunk)), 'frames');
        frames = cell(framesperchunk, 1);
    else
        frames{mod(i, framesperchunk)} = framedata;
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


end % unpack_wnmovie