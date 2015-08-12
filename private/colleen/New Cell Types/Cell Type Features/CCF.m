clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_opts.date='2009-04-13-7/'; % one slash at the end
run_opts.concatname='data000-cr'; % Name (or modified name) of run, no slashes

% Sometimes the data has two versions of the concate name
% run_opts.file_name = [run_opts.date, '/', 'data000-nwpca', '/',  'data000/data000'];
run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  run_opts.concatname];

% Full path to movie xml
run_opts.movie_spec='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-64x48.xml';


run_opts.save_location_root = '/Users/colleen/Desktop/Large Cell CCF/';
% Number of frames to use for generator signal as well as number of frames
% of the timecourse to display
run_opts.num_frames = 19;

% Number of bins to use for the nonlinearity graph
run_opts.num_bins = 10;

% How much padding to use for zooming in on STAs
params.padding = 7;

% Cell specification can be one cell type or multiple in a cell array.
% Use the same spelling/capitalization as the vision params file
% OFF large-2 in vision = OFF large 2
cell_specification = {'ON large type 1', 'ON large type 2', 'ON large type 3'};

%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Where to save
run_opts.filepath= [run_opts.save_location_root, run_opts.date, '/', run_opts.concatname, '/'];

% Used for labeling plot
run_opts.cell_type = cell_specification;

% Load Data
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opts.file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', run_opts.file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', run_opts.file_name, '.sta'];
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',1);
opt.load_sta_params.save_rf = 1; % Necessary to get the sta to load correctly
datarun=load_data(datarun,opt);

output = cell(size(cell_specification,2),1);
% Legend= cell(size(cell_specification,2),1);
fig = figure;
set(0,'DefaultAxesFontSize',14)
%         cmap = distinguishable_colors(size(cell_specification,2));
cmap = hsv(size(cell_specification,2));
count_leg = 1;
for type = 1:size(cell_specification,2)
    %     To figure out how many cells are in the cell type
    [cell_indices, cell_type_name] = get_cell_indices(datarun, cell_specification{type});
    cell_ids=datarun.cell_ids(cell_indices);
    
    
    N = length(cell_ids);
    if N > 1
        Legend{count_leg} = cell_specification{type};
        count_leg = count_leg+ 1;
        
        pairs_index = nchoosek(cell_indices,2);
        pairs_vision = nchoosek(cell_ids,2);
        
        output{type} = nan(size(pairs_vision));
        
        for i = 1:size(pairs_vision,1)
            % Scale location to pixels
            loc_cell1 = datarun.vision.sta_fits{pairs_index(i,1)}.mean* datarun.stimulus.stixel_width;
            loc_cell2 = datarun.vision.sta_fits{pairs_index(i,2)}.mean*datarun.stimulus.stixel_width;
            output{type}(i,1) = norm(loc_cell1 - loc_cell2);
            
            [time, ccf] = plot_ccf(datarun, [pairs_vision(i,1) pairs_vision(i,2)]);
            reduced_time = time(time >= -0.005 & time <= 0.005);
            reduced_ccf = ccf(time >= -0.005 & time <= 0.005);
            
            % figure; plot(time, ccf)
            [time_mean, ccf_mean] = plot_ccf(datarun, [pairs_vision(i,1) pairs_vision(i,2)], 'shuffle', 'stim', 'trial', 10);
            scale = mean(ccf_mean);
            output{type}(i,2) = trapz(reduced_ccf- scale)*mean(diff(reduced_time));
            
            
            H{type} =  plot(output{type}(:,1), output{type}(:,2), 'o', 'color', cmap(type,:));
            
            hold on
            
            
        end
        
    end
end

% Select the line axes properites of the first cell of each
% type
count = 1;
for j = 1:size(H,2)
    if ~isempty(H{j})
        inter(count) = H{j}(1);
        count = count+1;
    end
    
end

% Set the line properties and string properites of the legend
legend_handle = legend([inter] , Legend, 'Location', 'Northeast');


xlabel('Distance between cells (pixels)')
ylabel('Integrated cross correlation above chance')
%     title('Cross Correlation Characteristics among Cell Types')

spaces = strfind(cell_specification{2}, ' ');
title_phrase = cell_specification{2}(1:(spaces(2)-1));
suptitle({[run_opts.date, run_opts.concatname]; title_phrase})


% Make the directory to save in if it doesn't exist
if ~exist(run_opts.filepath)
    mkdir(run_opts.filepath)
end

% Save the figure to pdf
export_fig([run_opts.filepath, title_phrase], '-pdf')



