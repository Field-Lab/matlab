function h = stack_reslicer(stack, varargin)
% STACK_RESLICER    GUI to pick points in an image stack and then reslice the stack
% usage: handle = stack_reslicer(stack, opts)
%
% See STACK_POINT_PICKER and RESLICE_STACK for args.
%
% Imagine you have an image stack, but it's tilted or badly warped so that
% the RGC layer is not in any single section.  At the same time the max
% projection includes too much other stuff to be useful.  Run
% stack_reslicer and go through the sections, picking XY points for each 
% section where the RGCs are well represented in that section.  Then hit 
% reslice stack (from the pull down menu) to generate a new XY image that 
% picks the image points from the stack at the appropriate section for each 
% XY coordinate.
%
% See RESLICE_STACK for details.  This can take a long time for big stacks.
%
% 2010-05 phli
%

% TODO: 
%   Add a "save resliced" menu option
%   Grey interface and change pointer to spinner while reslicing
%

opts = inputParser();
opts.addParamValue('method', 'nearest');
opts.addParamValue('thickness', 1);
opts.addParamValue('verbose', true);
opts.addParamValue('gui', true);
opts.addParamValue('mesh_points', []);
opts.KeepUnmatched = true;
opts.parse(varargin{:});
unmatched = opts.Unmatched;
opts = opts.Results;

if ~isfield(unmatched, 'points') && isfield(stack, 'reslice_points')
    unmatched.points = stack.reslice_points;
end

h = stack_point_picker(stack, unmatched);

fig = getfig(h);
setappdata(fig, 'reslice_opts', opts);
setappdata(fig, 'mesh_points', opts.mesh_points);

% Add to API object
api = getappdata(fig, 'api');
api.gui_reslice_stack = @gui_reslice_stack;
api.show_resliced     = @show_resliced;
setappdata(fig, 'api', api);

% Add menu to reslice based on points
gh = guihandles(fig);
smenu = gh.smenu;
uimenu(smenu, 'Label', 'Relice stack',               'Tag', 'rsmenu',   'Callback', @gui_reslice_stack, 'Separator', 'on');
uimenu(smenu, 'Label', 'Show resliced',              'Tag', 'srmenu',   'Callback', @show_resliced);
ermenu = gui_data_export(fig, 'resliced', smenu);

if nargout < 1, clear h; end



function gui_reslice_stack(handle, ev) %#ok<INUSD>
fig = getfig(handle);
opts   = getappdata(fig, 'reslice_opts');
stack  = getappdata(fig, 'data');
points = getappdata(fig, 'points');

[XI, YI, ZI] = interp_scattered_stack_points(raw_stack(stack), points, opts);
mesh_points = {XI, YI, ZI};

resliced = reslice_stack(stack, mesh_points, opts);
setappdata(fig, 'resliced', resliced);
show_resliced(handle, []);

function show_resliced(handle, ev) %#ok<INUSD>
fig = getfig(handle);
resliced = getappdata(fig, 'resliced');
if ~isempty(resliced)
    figure;
    imshow(resliced);
end