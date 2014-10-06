function plot_pie_from_cone(datarun, cell_id, cone_location, varargin)


p = inputParser;

p.addParamValue('fig_or_axes', -1);
p.addParamValue('clear', true, @islogical);
p.parse(varargin{:})

params = p.Results


plot_axes = set_up_fig_or_axes(params.fig_or_axes, params.clear);
axes(plot_axes)

colormap_matrix = [1 0 0; 0 1 0; 0 0 1];


% plot time courses of stixels of interest
cell_index = get_cell_indices(datarun, cell_id);
datarun = load_sta(datarun, 'load_sta', cell_id, 'verbose', false);
datarun = get_sta_summaries(datarun, cell_id);
sta_size = size(datarun.stas.stas{cell_index});

% get a time course template and normalize the amplitude to green peak
tc_template = datarun.stas.time_courses{cell_index};
tc_template = mean(tc_template, 2);
norm_tc_template = tc_template ./ abs(min(tc_template));

% get tc of stixel and project along template
temp_sig_stix = false(sta_size(1:2));
temp_sig_stix(cone_location(1), cone_location(2)) = true;
temp_tc = time_course_from_sta(datarun.stas.stas{cell_index}, temp_sig_stix);
tc_projection = temp_tc' * norm_tc_template;

% make pie graph for m_one
pie(tc_projection * 100) % multiply by 100 to make sure sum(tc_projections) > 1
colormap(colormap_matrix)

