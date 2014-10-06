function sta = sta_from_java_stas(java_stas, cell_id, frames, independent)
% sta_from_java_sta     load and reshape one rrs sta from a java sta object
%
% usage:  sta = sta_from_java_sta(java_stas, cell_id, frames)
%
% arguments:   java_sta - java object generated by edu.ucsc.neurobiology.vision.io.STAFile(<sta_file_path>)
%               cell_id - cell id of the sta to get
%                frames - vector of which frames to get (choose from 0 through n-1)
%           independent - 't' for RGB sta, 'nil' for BW sta
%
% outputs:          sta - 4-d matrix of the STA
%
%
% 2008-10 gauthier
%
%



% VERIFY ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addRequired('java_stas',@isjava);
p.addRequired('cell_id',@isnumeric);
p.addRequired('frames', @(x)(any(strcmp(x,{':'})) || isnumeric(x)));
p.addRequired('independent', @(x)any(strcmpi(x,{'t','nil'})) );

% resolve user input and default values
p.parse(java_stas, cell_id, frames, independent);

% get params struct
params = p.Results;



% be sure cell id exists
cell_ids = java_stas.getIDList;
if ~isempty(setdiff(cell_id,cell_ids))
    error('cell id %d could not be found.',cell_id)
end



% get single-STA java object from the larger java object
java_sta=java_stas.getSTA(cell_id);


% extract STA data and convert to standard format

switch 1
    case 1
        sta = sta_from_java_sta(java_sta,'frames',params.frames,'independent',independent);
        
    case 2


        % format of rrs STAs (as pieced together by gauthier, 2008-03)
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
        height = java_stas.getHeight;
        width = java_stas.getWidth;

        % get list of desired frames
        frames = parse_frame_spec(frames, java_stas.getSTADepth);

        % set up variable to store data temporarily
        sta_temp = single(zeros(3,width,height,length(frames)));

        % load each desired frame
        for ff = 1:length(frames)
            % frame indexing begins at 0
            frame = java_sta.getFrame(frames(ff)-1).getBuffer;
            sta_temp(:,:,:,ff)=reshape(frame,3,width,height);
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
if strcmpi(independent,'nil')
    sta = sta(:,:,1,:);
end
