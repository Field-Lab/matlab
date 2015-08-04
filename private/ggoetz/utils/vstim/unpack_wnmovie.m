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

% Check if we're suppose to be unpacking a color or greyscale movie
bwmovie = strcmp(wnmovie.colorType().toString(), 'DEPENDENT');
if bwmovie
    nframes = nframes/3;
    if mod(framesperchunk, 3) ~= 0
        error('For BW movies, the number of frames per chunk should be a multiple of 3.')
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
    bufferdata = double(reshape(wnmovie.getFrame(i-1).getBuffer(), ...
        3, wnmovie.getWidth(), wnmovie.getHeight()));
    framedata = cat(3, squeeze(bufferdata(1,:,:)).', ...
        squeeze(bufferdata(2,:,:)).', squeeze(bufferdata(3,:,:)).');
    
    if bwmovie
        if mod(i, framesperchunk/3) == 0
            frames{framesperchunk-2} = squeeze(framedata(:,:,1));
            frames{framesperchunk-1} = squeeze(framedata(:,:,2));
            frames{framesperchunk} = squeeze(framedata(:,:,3));
            save(fullfile(outputpath, ...
                sprintf(filename, round(i*3/framesperchunk))), 'frames');
            frames = cell(framesperchunk, 1);
        else
            frames{(mod(i, framesperchunk/3)-1)*3 + 1} = squeeze(framedata(:,:,1));
            frames{(mod(i, framesperchunk/3)-1)*3 + 2} = squeeze(framedata(:,:,2));
            frames{(mod(i, framesperchunk/3)-1)*3 + 3} = squeeze(framedata(:,:,3));
        end
    else
        if mod(i, framesperchunk) == 0
            frames{framesperchunk} = framedata;
            save(fullfile(outputpath, ...
                sprintf(filename, i/framesperchunk)), 'frames');
            frames = cell(framesperchunk, 1);
        else
            frames{mod(i, framesperchunk)} = framedata;
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

end % unpack_wnmovie