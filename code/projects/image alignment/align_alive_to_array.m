function [array_tform, clicked_corners, ai] = align_alive_to_array(im, datarun_or_array_id, varargin)
% ALIGN_ALIVE_TO_ARRAY    Quick alignment method by clicking corners
% usage: [array_tform, clicked_corners, ai] = align_alive_to_array(im, datarun_or_arrayid, opts)
%
% inputs: datarun   
%         im        Image including electrodes to align to array coordinates
%
% opts: rig    []    If datarun does not specify datarun.piece.rig; needed if
%                    show_corners is true.
%
%       show_corners    true    Whether to help user by showing a plot of
%                               which order to click electrodes.
%
%       plot_radius     200     How many pixels to include around the
%                               initial rough alignment point when clicking
%                               the fine alignment point.
%
%       fig             []      Figure to plot in; if empty, uses gcf
%
% 2010-05 phli, ported from Jeff's 2009 script
% 2010-08 phli, removed hard dependency on datarun
%

% ToDo - waitforbuttonpress is not so nice in that it is triggered even if
% the clicks are on uibuttons or other figures.  Update to use normal
% callbacks.

opts = inputParser;
opts.addParamValue('rig', []);
opts.addParamValue('show_corners', true);
opts.addParamValue('plot_radius', 200);
opts.addParamValue('fig', []);
opts.parse(varargin{:});
opts = opts.Results;

% Plot image
fig = opts.fig;
if isempty(fig); fig = gcf; end
figure(fig); clf; imshow(im)

if isstruct(datarun_or_array_id)
    [datarun, ai] = load_array_info(datarun_or_array_id);
else
    datarun = struct();
    ai = load_array_info(datarun_or_array_id);
end

% Help user know which order to click?
if opts.show_corners
    if isempty(opts.rig)
        if isfield(datarun, 'piece') && isfield(datarun.piece, 'rig') && ~isempty(datarun.piece.rig)
            opts.rig = datarun.piece.rig;
        else
            error('Rig is not specified.  Try load_index, or set datarun.piece.rig, or pass rig as an option.');
        end
    end
    f2 = figure;
    show_alive_corner_orientation(ai.array_id, opts.rig);
end


% Pick corner points once roughly
fprintf('Click on corner electrodes in the following order: %s.\n(These clicks only need to be approximately correct.)\n', num2str(ai.corner_electrodes));
num_points = length(ai.corner_electrodes);
rough_points = zeros(num_points,2);
for pp=1:num_points
    waitforbuttonpress
    fum = get(gca,'CurrentPoint');
    rough_points(pp,:) = fum(1,1:2); % button down detected
end

% Do it again, this time zoomed in on the electrode
fprintf('click on corner electrodes again (this time the window will be zoomed in on the electrodes)\n')
clicked_corners = zeros(num_points,2);
for pp=1:num_points
    % plot on zoomed in region
    figure(fig);clf;imagesc(im);axis image;
    set(gca,'xlim', opts.plot_radius * [-1 1] + rough_points(pp,1),...
            'ylim', opts.plot_radius * [-1 1] + rough_points(pp,2))
    waitforbuttonpress
    fum = get(gca,'CurrentPoint');
    clicked_corners(pp,:) = fum(1,1:2); % button down detected
end
close(fig);
if ishandle(f2)
    close(f2);
end


% make array transformation
corner_electrode_positions = ai.positions(ai.corner_electrodes,:);
array_tform = cp2tform(clicked_corners, corner_electrode_positions, 'projective');