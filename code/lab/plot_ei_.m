function elec_handles = plot_ei_(ei, position, frame_number, varargin)
% PLOT_EI_  Plot a frame of an electrophysiological image
%
% frame_number==0 max(abs(ei))
%
% shlens@salk.edu 2006-07-17
% greschner new framework,
%   uses raw ei,
%   take full EI as input, if no framenumber is given or it's 0 -> produces max EI figure
%   shows positive and negative deflections in spike wave-form
%   no new figure
% 2010-02 phli
%   allows specifying specific RGB triple for each electrode
%


% SET UP OPTIONAL ARGUMENTS
p = inputParser;
p.addParamValue('plot_axes',gca);%
p.addParamValue('cutoff', 0.03);%
p.addParamValue('absolute_cutoff', true);%
p.addParamValue('max_scale', 1);% diameter of disks in units of median nearest neighbor spacing
p.addParamValue('scale', 3);%
p.addParamValue('elec_colors', []); % Give every electrode its own RGB triple
p.addParamValue('neg_color', [0 0 0]);%
p.addParamValue('pos_color', [0 0 0]);%
p.addParamValue('highlight', []);%electrode
p.addParamValue('highlight_color', [0 0 1]);%
p.addParamValue('boundary_scale', 1.1);%
p.addParamValue('alpha', true);%
p.addParamValue('pretty_axes', true);%
p.addParamValue('rotation', 0);%in deg
p.addParamValue('fliplr', 0);%
p.addParamValue('flipud', 0);%
p.addParamValue('zoom', 1);%
p.addParamValue('shift', [0 0]);%
p.addParamValue('label', false);%
p.addParamValue('labels', {}); % If given, use these as labels, otherwise just label by number
p.addParamValue('disk_points', 125);%
p.addParamValue('elec_callback', []);%
p.addParamValue('elec_spacing', infer_electrode_spacing(position));
p.addParamValue('axon', []);
p.addParamValue('axon_color', 'g');
p.addParamValue('sub_scale', 0);
p.addParamValue('stack_outline', []);
p.addParamValue('negative', 0);
p.addParamValue('log', 0);
p.addParamValue('start_frame', 1);
p.addParamValue('end_frame', []);
p.addParamValue('zlevel', []); % Helpful for plotting in front of / behind other data
% parse inputs
p.parse(varargin{:});
params = p.Results;


% if EI is empty or contains NaN, abort plotting
if isempty(ei) || any(isnan(ei(:)))
    text(mean(xlim(params.plot_axes)), mean(ylim(params.plot_axes)), 'invalid EI', 'HorizontalAlignment', 'center');
    return
end

% if EI is all 0, don't plot anything
if all(reshape(ei,[],1)==0)
    return
end

if isempty(params.end_frame)
    params.end_frame=size(ei,2);
end

%adjust position
if params.fliplr || params.flipud || params.rotation
    for i=1:length(position)
        if params.rotation
            [t1 t2]=cart2pol(position(i,1),position(i,2));
            t1=t1+(params.rotation/(180/pi));
            [position(i,1),position(i,2)]=pol2cart(t1, t2);
        end
        if params.fliplr
            position(i,1)=-position(i,1);
        end
        if params.flipud
            position(i,2)=-position(i,2);
        end
    end
end

position=position*params.zoom;
params.max_scale=params.max_scale*params.zoom;

position(:,1)=position(:,1)+params.shift(1);
position(:,2)=position(:,2)+params.shift(2);



% convert disk size from units of median nearest neighbor spacing to absolute units
params.max_scale = params.max_scale * params.elec_spacing/2;



%normalize
if ~params.negative
    ei = ei/max(abs(ei(:)));
else
    ei = ei/abs(min(ei(:)));    
end

% zero electrodes below cutoff
if params.absolute_cutoff
    % in the case of an absolute cutoff, use the cutoff value
    ei(abs(ei)<params.cutoff) = 0;
else
    % otherwise, make it a relative cutoff, and multiply by the scale first
    if ~params.negative
    	ei(abs(ei)*params.scale<params.cutoff) = 0;
    else
        ei(min(ei)*params.scale>params.cutoff) = 0;   
    end
    
end

% reduce the EI to a single frame

% if a nonzero frame number was specified
if exist('frame_number','var') && ~isempty(frame_number) && frame_number
    % use it
    ei_frame=ei(:,frame_number);
else
    % otherwise, set the value of each electrode to be the largest amplitude it ever attains, preserving the sign
    ei_frame = get_ei_max_frame(ei, params.start_frame, params.end_frame, params.negative);
end
params.scale;

% circle samples
circle_samples = [0:(2*pi/params.disk_points):2*pi 2*pi];
x_circle = cos(circle_samples);
y_circle = sin(circle_samples);


% draw a disk for each electrode not equal to 0
hold on;
for e = find(ei_frame ~= 0)'
    % Get color
    if ~isempty(params.highlight) && ismember(e, params.highlight)
        color = params.highlight_color;
    elseif ~isempty(params.elec_colors)
        color = params.elec_colors(e,:);
    else
        if (ei_frame(e) < 0)
            color = params.neg_color;
        else
            color = params.pos_color;
        end
    end

    % size of circle
    sd = params.max_scale * params.scale * (abs(ei_frame(e))-params.sub_scale);
    if sd > params.max_scale,
        sd = params.max_scale;
    end

    % stretch circle
    %S = [sd 0; 0 sd];
    %X = S * [x_circle; y_circle];
    X = sd * [x_circle; y_circle];
    
    % draw patch, or plot outline
    if params.alpha == false
        elec_handle = plot(params.plot_axes, X(1,:) + position(e,1), X(2,:) + position(e,2), 'Color', color);
    else
        elec_handle = patch(X(1,:) + position(e,1), X(2,:) + position(e,2), color, 'FaceAlpha', double(params.alpha), 'EdgeColor', color, 'Parent', params.plot_axes);
    end
    
    if ~isempty(params.zlevel)
        zdata = params.zlevel * ones(size(get(elec_handle, 'XData')));
        set(elec_handle, 'ZData', zdata);
    end
    
    % Helpful for click callbacks
    setappdata(elec_handle, 'num', e);
    
    if ~isempty(params.elec_callback)
        set(elec_handle, 'ButtonDownFcn', params.elec_callback);
    end
    
    elec_handles(e) = elec_handle;
end

% axis control
if params.pretty_axes
    %  set(gca,'XTick',[],'YTick',[]);
    %set(gca,'XDir','reverse','YDir','reverse');%make it agree with java display
    box on;
    axis equal;
    range=[min(position(:,1)) max(position(:,1)) min(position(:,2)) max(position(:,2))]*params.boundary_scale;
    axis(range);
end

% plot electrode numbers, if desired
if params.label
    for ee=1:size(position,1)
        if isempty(params.labels)
            labeltext = num2str(ee);
        else
            labeltext = params.labels{ee};
        end
        text(position(ee,1),position(ee,2),labeltext,'Color', params.pos_color,'FontSize',10,...
            'HorizontalAlignment','Center','VerticalAlignment','Bottom');
    end
end

% add axon, if provided
if ~isempty(params.axon)
    points = traced_cell_points(params.axon(1:2,:),params.axon(2:end,:));
    plot(points(:,1), points(:,2), params.axon_color);
end

if ~isempty(params.stack_outline)
    plot(params.stack_outline(:,1), params.stack_outline(:,2));
end

if nargout < 1; clear elec_handles; end