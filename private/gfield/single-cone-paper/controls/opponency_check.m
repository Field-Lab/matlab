%%
% compares the opponency estimated with 1x1 versus LMS stimuli.

% plantain 1x1
master_path{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
master_path{1,2} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard';

% plantain - LMS
slave_path{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data005/data005/data005';

s_contrast_factor = 1;
array_type = 519;

cell_types = {3,4};

% load master data
datarun{2} = load_data(master_path{1,1});
datarun{2} = load_params(datarun{2},struct('verbose',1));  
datarun{2} = load_ei(datarun{2}, 'all', struct('array_type', array_type));
datarun{2} = load_sta(datarun{2},'load_sta',[]);
datarun{2} = set_polarities(datarun{2});
datarun{2} = import_single_cone_data(datarun{2}, master_path{1,2});    
datarun{2}.cones.mosaic = make_mosaic_struct(datarun{2}.cones.centers);

%datarun{2} = get_sta_summaries(datarun{2}, cell_types);

% load slave data
datarun{1} = load_data(slave_path{1,1});
datarun{1} = load_params(datarun{1},struct('verbose',1));  
datarun{1} = load_neurons(datarun{1});
datarun{1} = load_ei(datarun{1}, 'all',struct('array_type', array_type));
datarun{1} = load_sta(datarun{1},'load_sta',cell_types);
datarun{1} = set_polarities(datarun{1});

marks_params.thresh = 3.0;
datarun{1} = get_sta_summaries(datarun{1}, cell_types,'marks_params', marks_params);




save_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/';
map_ei_classification_txt(datarun{1}, datarun{2},'master_cell_type',cell_types, 'slave_cell_type', cell_types,...
                        'classification_path', save_path, 'corr_threshold', 0.99);

datarun{2}.names.map_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/ei_map-from-.txt';
datarun = load_map(datarun);





% get opponency from single cone data
[mosaic_weights, selection, extras] = select_cone_weights(datarun{2}, {3},...
                                                'thresh', 0.03,...
                                                'radius', [0 8], 'polarity', 0,...
                                                'contiguity', false,'scale', 3.0,...
                                                'remove_cones', 'S');   


connectivity = mosaic_weights .* selection; % keep weights continuous valued

% the new datarun has S cone excluded from all fields of datarun.cone 
new_datarun = extras.new_datarun;

[num_cones, num_RGCs] = size(connectivity);


% normalize cone weights
for RGC = 1:num_RGCs
    connectivity(:,RGC) = connectivity(:,RGC) ./ sum(connectivity(:,RGC));
end

% compute purity of cone mosaic
purity_indices = compute_opponency_index(connectivity, new_datarun.cones.types);
figure(1)
hist(purity_indices,[-1:0.2:1])
axis([-1 1 0 20])
%%


% compute and store significant stixels
% thresh_params.select = 'thresh';
% thresh_params.thresh = 2.0;
% datarun{1} = get_significant_stixels(datarun{1}, cell_types);
% 
% 
% % compute rfs from significant stixels and stas
% cell_indices = get_cell_indices(datarun{1}, cell_types);
% for rgc = 1:length(cell_indices)
%     rf_params.frames = ':';
%     rf_params.sig_stixels = datarun{1}.stas.marks{cell_indices(rgc)};
%     rf = rf_from_sta(datarun{1}.stas.stas{cell_indices(rgc)}, rf_params);
%     datarun{1}.stas.rfs{cell_indices(rgc)} = rf;
% end


OIs = plot_dkl_cone_weights(datarun{1}, {3}, 'method', 'area_around_peak', 'window_size', 2,...
                    'S_contrast_factor', s_contrast_factor, 'foa', 5, 'sigma_cutoff', 3);
figure(2)
hist(OIs, [-1:0.2:1])
axis([-1 1 0 20])

%%

rand_indices = randperm(length(OIs));

figure(3)
plot(OIs, purity_indices, 'ko')
hold on
plot([-1 1], [-1 1], 'k')
axis square
xlabel('LMS')
ylabel('1x1')

corrcoef(OIs, purity_indices)
axis([-1 1 -1 1])

print(3, '~/Desktop/LMS-comp.pdf', '-dpdf')












