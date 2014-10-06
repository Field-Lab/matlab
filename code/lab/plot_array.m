function h = plot_array(datarun,color_spec,T,label)
% plot_array     plot the outline of the array
%
% usage:  h = plot_array(datarun,color_spec,T)
%
% arguments:     datarun - datarun struct with fields
%                           datarun.piece.corners
%            color_spec - color to plot with, if not specified default is black
%                     T - transformation to apply to the corners before plotting
%
%
% 2010-03  gauthier
%


if ~exist('color_spec','var')
    color_spec = 'k';
end

if isfield(datarun,'piece') && isfield(datarun.piece,'corners')
    corners = datarun.piece.corners;
    if exist('T','var') && ~isempty(T)
        corners = tformfwd(T,corners);
    end
    h = plot(corners(:,1),corners(:,2),'color',color_spec);
else
    h = [];
    warning('could not plot array outline because datarun.piece.corners does not exist.  see compute_monitor_to_array_transformation.');
end

if exist('label','var') && label
    if isfield(datarun.piece,'corner_electrodes')
        for cc = 1:length(datarun.piece.corner_electrodes)
            text(corners(cc,1),corners(cc,2),sprintf('%d',datarun.piece.corner_electrodes(cc)),...
                'Color',color_spec,'FontSize',15,'HorizontalAlignment','Center','VerticalAlignment','Bottom')
        end
    else
       warning('corner electrodes have not been identified.');
    end
end