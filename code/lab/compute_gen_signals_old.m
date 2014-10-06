function gen_signals = compute_gen_signals_old(strfs, movie, varargin)
% verbose       true                        show output
% rois          []                         	YxXxN matrix specifying which stixels to use for calculation
%                                               if empty, use all stixels for all cells
%                                               can be sparse and/or boolean
% start_stim    1                           first movie frame used in analysis
% end_stim      movie.size                  last movie frame
% num_frames    1                           how many temporal frames are in each STRF

opts = inputParser();
opts.addParamValue('verbose',       true );
opts.addParamValue('rois',          [] );
opts.addParamValue('start_stim',    1,          @(x)numel(x)==1 && x>=0);
opts.addParamValue('end_stim',      movie.size, @(x)numel(x)==1 && x>=0 && x <=movie.size);
opts.addParamValue('num_frames',    1);
opts.addParamValue('isRGB',         true,       @islogical)
opts.parse(varargin{:});
opts = opts.Results;

% note stimulus size
field_width = movie.getWidth;
field_height = movie.getHeight;

% identify numbers of things
start_stim = opts.start_stim;
end_stim = opts.end_stim;
num_frames = opts.num_frames;
num_stims = end_stim - start_stim + 1;
num_cells = length(strfs);

% first_stim is the earliest stimulus frame for which a generator signal value can be computed
first_stim = start_stim + num_frames - 1;


% GET GENERATOR SIGNAL VALUES

% show output
if opts.verbose
    fprintf('Computing generator signal values for frames %d to %d...',start_stim,end_stim)
    T=text_waitbar;
    start_time = clock; % note when it started
end


% initialize storage variables
gen_signals = zeros(num_stims-num_frames+1,num_cells);


% cycle through each stimulus
for ss = start_stim:end_stim

    if opts.verbose
        T=text_waitbar(T,(ss-start_stim)/num_stims - 0.01);
    end
    
    
    % get new frame
    STAFrame = movie.getFrame(ss-1);
    new_frame = permute(reshape(STAFrame.getBuffer,3,field_width,field_height),[3 2 1]) - .5;

    % if first_stim hasn't been reached yet
    if ss < first_stim
        % load this frame
        stims(:,ss-start_stim+1+1) = reshape(new_frame,[],1);
        % keep going
        continue
    end
        
    % shift current frames over
    stims(:,1:end-1) = stims(:,2:end);
    % assign new frame to the final position
    stims(:,end) = reshape(new_frame,[],1);


    % get generator signal values for each cell
    for cc = 1:num_cells
        % get STRF
        strf = strfs{cc};

        % get stimulus in relevant region
        if opts.isRGB
            roi = repmat(full(opts.rois(:,cc)),3,1);
        else
            roi = repmat(full(opts.rois(:,cc)),1,1);
        end
        stim = reshape(stims(roi,:),[],1);

        % take dot product
        stim_num = ss - first_stim + 1;
        gen_signals(stim_num,cc) = stim'*strf;
    end


end


% display how long it took
if opts.verbose
    fprintf('done (%0.1f seconds)\n',etime(clock,start_time));
end