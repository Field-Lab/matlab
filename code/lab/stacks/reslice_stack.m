function resliced = reslice_stack(stack, points, varargin)
% RESLICE_STACK    Generate a new XY slice through an image stack with flexibly picked Z locations
% usage: resliced = reslice_stack(stack, points, opts)
%
% Stack can be either a 3D monochrome XYZ matrix, or a 4D X-Y-RGB-Z.
%
% Points should be in the format output by stack_point_picker: a cell array
% holding an Nx2 matrix of point XY coordinates for each section.
%
% opts  method      'nearest'   Set interpolation method; either 'nearest' or 'linear'.
%
%       channels    []          Set specific channel numbers to use.  Throw
%                               away useless channels without doing time
%                               consuming calculations on them.
%
%       verbose     true
%
%
% Imagine you have an image stack, but it's tilted or badly warped so that
% the RGC layer is not in any single section.  At the same time the max
% projection includes too much extra stuff to be useful.  RESLICE_STACK
% allows you to specify a series of XYZ points that indicate which stack 
% section has the desired information for each XY coordinate.
%
% The XYZ points are specified in a kind of funny format, to interact
% directly with the STACK_RESLICER GUI (should probably add another
% abstraction layer in here).
%
% The scattered XYZ points marking desired areas of the stack are first
% interpolated into a complete XY meshgrid giving the Z coordinate for eacy
% XY coordinate.  Then the stack images are interpolated to those Z values
% to generate the resliced XY image.
%
%
% 2010-05 phli
%

opts = inputParser;
opts.addParamValue('method', 'nearest');
%opts.addParamValue('thickness', 1);
opts.addParamValue('channels', []);
opts.addParamValue('verbose', true);
opts.addParamValue('progress', true);
opts.addParamValue('batchsize', 500);
opts.addParamValue('gui', false);
opts.KeepUnmatched = true;
opts.parse(varargin{:});
opts = opts.Results;

% For now, no real gui feedback so just turn on stdout
if opts.gui; opts.verbose = true; end


% Get formats
stack_class = class(stack);
if strcmp(stack_class, 'struct')
    stack = raw_stack(stack);
    stack_class = class(stack);
end


% Determine whether points passed are in meshed XI, YI, ZI format or still 
% in the scattered point format from stack_point_picker
psize = size(points{1});
if ndims(psize) == 2 && psize(1) == size(stack,1) && psize(2) == size(stack,2)
    XI = points{1};
    YI = points{2};
    ZI = points{3};
else
    [XI, YI, ZI] = interp_scattered_stack_points(stack, points, opts);
end


% Calculate batch processing bounds
xstart = 1:opts.batchsize:size(stack,1);
xend = opts.batchsize:opts.batchsize:size(stack,1);
ystart = 1:opts.batchsize:size(stack,2);
yend = opts.batchsize:opts.batchsize:size(stack,2);
if length(xend) < length(xstart), xend(end+1) = size(stack,1); end
if length(yend) < length(ystart), yend(end+1) = size(stack,2); end


% Process in batches to avoid memory swapping
if opts.verbose, disp('Reslicing...'); end
if opts.progress, for j = 1:length(xstart), for k = 1:length(ystart), fprintf('='); end; end; fprintf('\n'); end
resliced = zeros([size(stack, 1), size(stack, 2), size(stack, 3)], stack_class);
for j = 1:length(xstart)
    xb = xstart(j):xend(j);
    
    for k = 1:length(ystart)
        yb = ystart(k):yend(k);
        xi = XI(1:length(xb),1:length(yb));
        yi = YI(1:length(xb),1:length(yb));

        % 3D is mono, 4D is separate RGB color channels
        if ndims(stack) == 3
            resliced(xb,yb) = cast(reslice_mono_stack(stack(xb,yb,:), xi, yi, ZI(xb,yb), opts, 'verbose', false), stack_class);
        elseif ndims(stack) == 4
            
            % Do each color channel independently
            if isempty(opts.channels), opts.channels = 1:size(stack, 3); end
            for i = opts.channels
                stack_channel = squeeze(stack(xb,yb,i,:));
                if ~any(stack_channel(:)), continue; end
                resliced(xb,yb,i) = cast(reslice_mono_stack(stack_channel, xi, yi, ZI(xb,yb), opts, 'verbose', false), stack_class);
            end
        else
            error('RESLICE_STACK:argerr', 'Invalid stack argument');
        end
        if opts.progress, fprintf('.'); end
    end
end
if opts.progress, fprintf('\n'); end