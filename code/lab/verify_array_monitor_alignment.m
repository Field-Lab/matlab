function verify_array_monitor_alignment(datarun,varargin)
% verify_array_monitor_alignment     Verify alignment points are well chosen
%
% usage:  verify_array_monitor_alignment(datarun, <params>)
%
% arguments:     datarun - datarun struct, it is assumed that the function
%                           compute_monitor_to_array_transformation was already run
%            varargin - struct or list of optional parameters (see below)
%
%
% optional params, their default values, and what they specify:
%
% which_pm        	'pm32'    	which image to use for alignment (for photographic mapping only)
% figure        	0          	figure or axes to plot in. if 0, make new figure.
%
%
% 2010-02  gauthier
%



% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% note valid names for alignment images
valid_image_names = load_alignment_image;

% specify list of optional parameters
p.addParamValue('which_pm','pm32',@(x)any(strcmpi(x,valid_image_names)));
p.addParamValue('figure',0);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



if ~isfield(datarun,'piece') || ~isfield(datarun.piece,'T_monitor_to_array')
    error('Monitor alignment was not performed, because datarun.piece.T_monitor_to_array does not exist.')
end

% figure out which kind of mapping was done
if isfield(datarun.piece,'photographic_mapping')
    mapping_type = 'photographic';
elseif isfield(datarun.piece,'array')
    mapping_type = 'clicked';
else
    error('It appears monitor alignment was not performed.')
end



switch mapping_type
    case 'photographic'
        % PLOT VERIFICATION OF PHOTOGRAPHIC MAPPING


        % GET STUFF TO PLOT

        % get photographic mapping image
        %pm_image = datarun.piece.photographic_mapping.images.(params.which_pm);
        pm_image = imread([server_data_path '/' datarun.piece.(sprintf('photographic_mapping_%s',params.which_pm))]);

        % get alignment image
        ai = load_alignment_image('pm',params.which_pm);

        % get electrode positions
        array_info = load_array_info(datarun,2);
        % transform electrodes to base coordinates
        epb=tforminv(datarun.piece.photographic_mapping.T_base_to_array,array_info.positions);



        % PLOT PM IMAGE

        % get pm in array_coordinates
        [ai_array.im,ai_array.xdata,ai_array.ydata] = ...
            imtransform(ai.im,fliptform(datarun.piece.photographic_mapping.T_base_to_monitor),...
            'udata',ai.xdata,'vdata',ai.ydata,'xyscale',1,...
            'xdata',[1 size(pm_image,2)],'ydata',[1 size(pm_image,1)] );

        % set up figure
        set_up_fig_or_axes(params.figure);
        fig_num = gcf;

        % plot it
        figure(fig_num);
        subplot(131)
        imagesc(ai_array.im,'xdata',ai_array.xdata,'ydata',ai_array.ydata)
        axis image; colormap gray;hold on
        title('predicted photographic mapping image')

        % add electrodes
        plot(epb(:,1),epb(:,2),'.')



        % PLOT PHOTOGRAPH

        figure(fig_num);

        subplot(132)
        imagesc(pm_image)
        axis image; hold on
        title('actual photographic mapping image')

        % plot electrodes
        plot(epb(:,1),epb(:,2),'.')

        % link the axes
        linkaxes(get(fig_num,'child'))
        
        
        
        % PLOT ARRAY

        % get photographic mapping image
        %pm_image = datarun.piece.photographic_mapping.images.array_edges;
        pm_image = imread([server_data_path '/' datarun.piece.photographic_mapping_array_edges]);

        % get pm in array_coordinates
        [ai_array.im,ai_array.xdata,ai_array.ydata] = ...
            imtransform(ai.im,fliptform(datarun.piece.photographic_mapping.T_base_to_monitor),...
            'udata',ai.xdata,'vdata',ai.ydata,'xyscale',1,...
            'xdata',[1 size(pm_image,2)],'ydata',[1 size(pm_image,1)] );


        figure(fig_num);

        subplot(133)
        imagesc(pm_image)
        axis image; hold on
        title('actual photographic mapping image')

        % plot electrodes
        plot(epb(:,1),epb(:,2),'+')

        % link the axes
        linkaxes(get(fig_num,'child'))


    case 'clicked'
        % PLOT VERIFICATION OF CLICKED POINT MAPPING

        % set up figure
        set_up_fig_or_axes(params.figure);
        fig_num = gcf;
        figure(fig_num);
        
        % plot background image for monitor
        mon_x = datarun.stimulus.monitor_x;
        mon_y = datarun.stimulus.monitor_y;
        image(64*ones(mon_y,mon_x));axis image; hold on;colormap gray
        
        % plot clicked points
        plot(datarun.piece.array(:,1),datarun.piece.array(:,2),'r.','Markersize',25)
        
        % plot interpolated map
        plot(datarun.piece.corners(:,1),datarun.piece.corners(:,2),'k.-')

        title(sprintf('red dots are clicked points, black dots and line are estiamted array position.\nplotted in monitor coordinates.'))
end

