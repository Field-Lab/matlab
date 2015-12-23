function average_frame = compute_average_white_noise_frame(source,varargin)
% compute_average_white_noise_frame     Just like the name says
%
% usage:  frame = compute_average_white_noise_frame(datarun,varargin)
%         frame = compute_average_white_noise_frame({'.../x.movie','.../x.globals'},varargin)
%
% arguments:   source - datarun struct with datarun.names.rrs_prefix
%                           or cell array of two strings: path to movie file, path to globals file
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     frame - average frame
%
%
% optional params, their default values, and what they specify:
%
% verbose           false               show output
%
%
% 2008-10 gauthier, machado
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of required parameters
p.addRequired('source',@(x)isstruct(x)|iscellstr(x));

% specify list of optional parameters
p.addParamValue('verbose', true);

% resolve user input and default values
p.parse(source,varargin{:});

% get params struct
params = p.Results;






%make average movie frame
%you need the latest version of vision on your javapath

if params.verbose;
    name = 'Computing Average White Noise Movie Frame';
    fprintf('%s ...\n', name);
    start_time = clock;
end

%get path to .movie and .globals
switch class(source)
    case 'struct' % datarun struct, use rrs prefix to getthe name
        movieName = [source.names.rrs_prefix '.movie'];
        globalsName = [source.names.rrs_prefix  '.globals'];

    case 'cell' % names are given in a cell array
        movieName = source{1};
        globalsName = source{2};
end

%make movie object
mf = edu.ucsc.neurobiology.vision.stimulus.WhiteNoiseMovie(movieName,globalsName);

%get number of frames
nFrames = mf.size;

%there are always three guns
nGuns = 3;

% create text progress bar
if params.verbose;
    T=text_waitbar('Processsing Frames');
end

%calculate size of each frame
w = mf.getWidth; h = mf.getHeight;
fr_reshaped = zeros(h,w,nGuns);

%get each frame
for frame=0:nFrames-1
    fr = mf.getFrame(frame).getBuffer;
    for gun=1:nGuns
        fr_reshaped(:,:,gun) = fr_reshaped(:,:,gun) + reshape(fr(gun:3:h*w*nGuns),w,h)';
    end
    T=text_waitbar(T,frame/nFrames);
end
fprintf('\n');

average_frame = fr_reshaped ./ nFrames;

%plot the frame
figure(1);clf;imagesc(frame); axis image;
fprintf('%s ... done (%0.1f seconds)\n', name, etime(clock,start_time));


% 
% %get each frame
% for frame=0:nFrames-1
%     if frame==0
%         % get size of a frame
% 
%         %calculate size of each frame
%         fSize = h*w*nGuns;
%         fr = zeros(fSize,1);
%     end
% 
%     % add to accumulating sum
%     fr = fr + mf.getFrame(frame).getBuffer;
%     
%     % show progress
%     if params.verbose;
%         T=text_waitbar(T,frame/nFrames);
%     end
% end
% fprintf('\n');
% 
% % divide by number of frames to average
% fr = fr/nFrames;
% 
% %reshape
% average_frame = zeros(h,w,nGuns);
% for gun=1:nGuns
%     average_frame(:,:,gun) = reshape(fr(gun:3:h*w*nGuns),h,w);
% end
% 
% if params.verbose;
%     fprintf('%s ... done (%0.1f seconds)\n', name, etime(clock,start_time));
% end
% 
