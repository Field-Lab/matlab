% compare hot spots to parasols 


clear
close all
% dbstop if error
dataparam.date='2016-04-21-1';
dataparam.concatname='data005';
dataparam.jitter_concatname='data006-cf/edited/data006-cf';
dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-2-0.48-22222-119.5.xml';
dataparam.stixel_size = 8;
dataparam.seed = 22222;
fitparam.num_frames = 30;
frame_width = 640/dataparam.stixel_size;
frame_height = 320/dataparam.stixel_size;
stixels_per_frame = frame_width*frame_height;

num_colors =3;
dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];

jitter_id = 6513;
load(['/Volumes/Lab/Users/crhoades/Jitter/', dataparam.date, '/', dataparam.jitter_concatname, '/Cell ', num2str(jitter_id), '.mat'])


% list specific cell (1), or run for a whole cell type (0)
select_cells = 1;
if select_cells == 1
   dataparam.cell_specification = [4 5538 5568 5641 5674 6032 6140 6183 6424 6496 6663]; %ON parasol

end
dataparam.cell_type = {'all'};
%% END OF INPUT


% Right Movie
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name_right, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name_right, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name_right, '.sta'];


%% Load Data2

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false, 'load_sta', 0, 'load_sta_params', 1, 'load_all',true);
opt.load_sta_params.save_rf = 1;
% opt.load_sta_params.frames = 1:fitparam.num_frames;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);
datarun = load_sta(datarun, struct('load_sta',dataparam.cell_specification));
scale = 640/datarun.stimulus.field_width; 
% scale = 1;

sta = permute(temp, [2 1 3 4]);
[junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,4)),1));

% normalize STA color
sta = norm_image(sta);
figure; 
image(sta(:,:,:,start_index));
axis image

hold on 


% plot_axes = plot_rf_summaries(datarun, dataparam.cell_specification)
[cell_indices] = get_cell_indices(datarun, dataparam.cell_specification);
  % figure
    for cell_index = 1:length(cell_indices)

        % get the fit
        the_fit = datarun.stas.fits{cell_indices(cell_index)};

        % skip if doesn't exist
        if isempty(the_fit);continue;end

        % get center
        ctr = the_fit.mean*scale;
        % get center radius
        rad = the_fit.sd*scale;

        % get points of an ellipse with these parameters

        [X,Y] = drawEllipse([ctr rad the_fit.angle]);
        
        % if any are NaN, skip
%         if any(isnan([X Y]));continue;end
        
        % transform to desired coordinates
%         [X, Y] = tformfwd(coord_tform, X, Y);

        % plot the points and fill
        
      
       plot(X,Y,'Color',[0 0 0 ], 'LineWidth', 1);
        hold on
    end
    set(gca, 'ydir', 'reverse')
    
    title({[dataparam.date, ' ', dataparam.jitter_concatname];['Cell ', num2str(jitter_id)]});
    axis off;
