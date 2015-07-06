function gen_signals = compute_gen_signals(strfs, movie, varargin)
% COMPUTE_GEN_SIGNALS   
% usage: gen_signals = compute_gen_signals(strfs, movie, opts)
%
% Optimized from the old version.  Closely related to the functions for
% Jeremy Freeman's cone generator signal calculations, CALC_CONE_INPUTS,
% etc.  Potentially should be merged.  Ideally this should be made abstract
% enough to handle multiple movie types; this could also follow the
% function handle / movie object approach of CALC_CONE_INPUTS.
%
% Still additional optimization possibilities here, but this covers the
% major issues.  The remaining large improvement possibility would be to
% make the Java movie capable of returning mono frames instead of RGB
% frames for monochromatic movies.  Probably worth testing.
%
% opts:
%   verbose       true          show output
%   start_stim    1             first movie frame used in analysis
%   end_stim      movie.size    last movie frame
%   num_frames    1             how many temporal frames are in each STRF
%   isRGB         true          is movie RGB  
%
% 2013-06 phli, optimized from compute_gen_signals_old
%
% See also: COMPUTE_GEN_SIGNALS_OLD, GET_SNLS, CALC_CONE_INPUTS
%

opts = inputParser();
opts.addParamValue('verbose',       true );
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
start_stim = opts.start_stim -1;
end_stim = opts.end_stim;
num_frames = opts.num_frames;
num_stims = end_stim - start_stim;
num_cells = length(strfs);


% Convert strf into visbuf format
for cc = 1:num_cells
    strfs{cc} = stavect2visbuf(strfs{cc}, movie.getWidth(), movie.getHeight(), opts.isRGB);
end


if opts.verbose
    fprintf('Computing generator signal values for frames %d to %d...', start_stim, end_stim)
    T = text_waitbar();
    start_time = clock; % note when it started
end

% Initialize
diag_gen_signals = cell(num_cells,1);
for cc = 1:num_cells
    diag_gen_signals{cc} = zeros(num_stims, num_frames);
end

% Step through movie frame by frame
% The gen_signal for each frame of the STRF is calculated separately for
% each movie frame here (diag_gen_signals).  The diag_gen_signals matrix is
% then summed along its diagonals below to calculate the total gen_signal
% for the whole STRF for each movie frame.
%for ss = start_stim:end_stim
for ss = 1:num_stims
    if opts.verbose, T = text_waitbar(T, ss / num_stims - 0.01); end
    
    stimbuf = movie.getFrame(ss-1+start_stim).getBuffer();
    if ~opts.isRGB, stimbuf = stimbuf(1:3:end); end

    % Get generator signal values for each cell
    for cc = 1:num_cells
        diag_gen_signals{cc}(ss,:) = double(stimbuf')*strfs{cc};
    end
end

if opts.verbose
    fprintf('done (%0.1f seconds)\n',etime(clock,start_time));
end


% Sum across strf frames 
% The result of the above is for each cell a series of vectors giving the
% generator signal across the whole movie for each separate frame of the
% STRF.  We want to sum across STRF frames for the total generator signal,
% which means summing along diagonals.  conv2 with an identity matrix is an
% efficient way to get Matlab to sum down all the diagonals.
gen_signals = zeros(num_stims - num_frames + 1, num_cells);
diag_sum_kernel = eye(num_frames);
% size(diag_gen_signals{cc});
% size(diag_sum_kernel)
for cc = 1:num_cells
    gen_signals(:,cc) = conv2(diag_gen_signals{cc}, diag_sum_kernel, 'valid');
end


% Unbias gen_signals
% The frames from Java movie are biased by +0.5.  This could be corrected
% in the main loop, but more efficient to compensate here.
for cc = 1:num_cells
    gen_signals(:,cc) = gen_signals(:,cc) - 0.5*sum(strfs{cc}(:));
end