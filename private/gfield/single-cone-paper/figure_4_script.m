
% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00';

% apricot
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2009-04-13-5/data005-s3600-s7200/data005/data005';
path_and_name{1,2} = '2009-04-13-5_data005-s3600-s7200_data005_data005-bayes-msf_15.00';


%'blueberry';
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';
path_and_name{1,2} = '2008-08-26-2_data001-s6369-s9551_data001-s6369-s9551-bayes-msf_15.00';


% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);


plot_dkl_cone_weights(datarun, {3,4}, 'method', 'cones')

print(1,'~/Desktop/figure1.pdf', '-dpdf')

% for plantain
% L center
% 813 ; nice ON cell
% 571 ; OK
% 766 ; nice

% M center
% 1099; nice off cell
% 976 ; nice
% 1666; nice

plot_dkl_cone_weights(datarun, 1099)

plot_cell_sampling(datarun, 1099, 'fig_or_axes', 10, 'line_width', [realmin 2.0],...
                    'thresh', 2.5, 'polarity', 0, 'contiguity', false, 'plot_radius', [0 6], 'cone_size', 15)
axis square
axis off

print(10, '~/Desktop/spider.pdf', '-dpdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot DLK cone weights from LMS isolating stimulus

% plantain 
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data005/data005/data005';
s_contrast_factor = 1;

% peach
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data004/data004/data004';
s_contrast_factor = 1;

% apricot
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2009-04-13-5/data006/data006';
s_contrast_factor = 4.5;

%blueberry
path_and_name{1,1} ='/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data006/data006/data006';
s_contrast_factor = 1;

%%%%%%%%%%%

cell_types = {3,4};

% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',cell_types);
datarun = set_polarities(datarun);
cell_indices = get_cell_indices(datarun, cell_types);

% compute and store significant stixels
thresh_params.select = 'thresh';
thresh_params.thresh = 2.0;
datarun = get_significant_stixels(datarun, cell_types);


% compute rfs from significant stixels and stas
for rgc = 1:length(cell_indices)
    rf_params.frames = ':';
    rf_params.sig_stixels = datarun.stas.marks{cell_indices(rgc)};
    rf = rf_from_sta(datarun.stas.stas{cell_indices(rgc)}, rf_params);
    datarun.stas.rfs{cell_indices(rgc)} = rf;
end


OIs = plot_dkl_cone_weights(datarun, cell_types, 'method', 'area_around_peak', 'window_size', 2,...
                    'S_contrast_factor', s_contrast_factor, 'foa', 1, 'sigma_cutoff', 3);

%title('plantain')
%print(1, '~/Desktop/cone_weights.pdf', '-dpdf')








