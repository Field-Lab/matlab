function sta = sta_from_java_sta(java_sta, varargin)
% sta_from_java_sta     convert java sta object to an STA in standard form
%
% usage:  sta = sta_from_java_sta(java_sta, varargin)
%
% arguments: java_sta - java sta object
%            varargin - struct or list of optional parameters (see below)
%
% outputs:        sta - 4-d matrix of the STA
%
%
% optional params, their default values, and what they specify:
%
% frames           	':'             which frames, see parse_frame_spec
% independent       't'             't' or 'nil'
%
%
%
% 2008-12 gauthier
% 2010-02 phli, changed to case 1; faster, cleaner
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('frames', ':');
p.addParamValue('independent', 't');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION


% the format of rrs STAs was pieced together by gauthier, 2008-03
% (this code accurately reconstructs the STAs, but it would be nice to paste
% some more detailed information about the STA format from the RRS documentation)
%
%
% each frame is stored as a single vector in the buffer
% the vector is a reshaped list of values from the original 3-d matrix
% sort order is y, x, and color
%
% each consecutive triple gives the r, g, and b values at one stixel
% for an rgb STA, these values will be different.
% for a bw STA, they will be identical.
%



% get STA size
height = java_sta.getHeight;
width = java_sta.getWidth;

% get list of desired frames
frames = parse_frame_spec(params.frames, java_sta.getSTADepth);

switch 1
    case 1
        sta = zeros(height, width, 3, length(frames), 'single');
        for i = 1:length(frames)
            sta(:,:,:,i) = sta_frame_from_java_sta_frame(java_sta.getFrame(frames(i)-1), height, width);
        end
    case 2
        % set up variable to store data temporarily
        sta_temp = single(zeros(3,width,height,length(frames)));
        
        
        % load each desired frame
        for ff = 1:length(frames)
            % frame indexing begins at 0
            frame = java_sta.getFrame(frames(ff)-1).getBuffer;
            sta_temp(:,:,:,ff) = reshape(frame,3,width,height);
        end
        
        clear java_sta  % does this clear it from memory?  it doesn't appear to have any method for clearing itself
        
        % initialize variable where STA will be stored
        sta = single(zeros(height,width,3,length(frames)));
        
        % go through each frame
        for ff = 1:size(sta_temp,4)
            % go through each color
            for rr = 1:3
                % reshape sta to put it into standard format
                sta(:,:,rr,ff) = flipud(rot90(squeeze(sta_temp(rr,:,:,ff))));
            end
        end
end

% if BW, remove duplicate copies in the color dimension
if strcmpi(params.independent,'nil')
    sta = sta(:,:,1,:);
end