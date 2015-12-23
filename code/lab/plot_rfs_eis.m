function plot_rfs_eis(datarun,cell_spec)
% plot_rfs_eis     plot many RFs with EIs in a slider window
%
% usage:  plot_rfs_eis(datarun,cell_spec)
%
% arguments:      datarun - datarun struct
%               cell_spec - 
%
%
% 2010-04  gauthier
%



% get cell ids
cell_ids = get_cell_ids(datarun,cell_spec);

% create slider control
figure
ha = make_loop_slider_list(1, 1, length(cell_ids), {@slider_plot, datarun,cell_ids});

% plot once before any clicks
slider_plot(ha, [], datarun,cell_ids);


function slider_plot(handle, event, datarun,cell_ids) %#ok<INUSL>
% display one frame of the STA

% get the slider position
cc = round(get(handle,'Value'));
cla;

% get cell id
cell_id = cell_ids(cc);

% plot RF in current axes
plot_rf(datarun,cell_id,'foa',-1,'array',true)

% add EI
hold on;
plot_ei(datarun,cell_id,'coordinates','sta','pretty_axes',0,'alpha',0,...
    'neg_color','b','pos_color','r','scale',1)
