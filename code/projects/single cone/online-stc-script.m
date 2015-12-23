% For a datarun with single cone RFs, compute the input to each cone.  
%
% set time_offset to change the time offset between spikes and stimuli.  Usually -1 is best.
%
% S = # stimuli
% C = # cones
% N = # cell ids
%
% Results are stored in three variables:
%
% cone_inputs - S x C matrix, each entry a cone activation value
% spike_times - N x S matrix, each entry number of spikes during that stimulus
% spike_rate  - length N cell, each entry a list of spike times (in units of stimulus number)
%
% cone_inputs are saved to a text file located at "save_path"
%

%data_path = '/snle/acquisition/2010-10-18-1';
data_path = '/snle/lab/Experiments/Array/Analysis/2010-09-24-1/streamed/data006/data006';
movie_path = '/snle/acquisition/movie-xml2/BW-1-6-0.48-11111-320x320.xml';
nickname = 'orange-6';

% SET PARAMETERS
start_time = 0;
end_time = 1799;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATARUN
datarun = load_data(data_path);

% load spikes times and trigger times
datarun = load_neurons(datarun);

% load java object of movie
datarun = load_java_movie(datarun, movie_path); 

datarun.names.nickname = nickname;

% load cone weights matrix
load([single_cone_path datarun.names.nickname '/Wc.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% refresh time
refresh_time = datarun.stimulus.java_movie.getRefreshTime/1000;

% compute start and stop times in terms of stimulus number
start_stim = floor(1+start_time/(datarun.stimulus.java_movie.getRefreshTime/1000));
end_stim = floor(1+end_time/(datarun.stimulus.java_movie.getRefreshTime/1000));
stims = start_stim:end_stim;

% note stimulus size
field_width = datarun.stimulus.java_movie.getWidth;
field_height = datarun.stimulus.java_movie.getHeight;


datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = import_single_cone_data(datarun, datarun.names.nickname);
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);
datarun = get_sta_summaries(datarun, {1,2}, 'keep_stas', false);


%%
% IDENTIFY SPIKES IN STIMULUS BINS


% get list of cells
cell_indices = get_cell_indices(datarun,'all');

% initialize storage variables
spike_rate = zeros(length(cell_indices),length(stims));
spike_times = cell(length(cell_indices),1);


% go through each cell
for cc = 1:length(cell_indices)

    cell_index = cell_indices(cc);

    % bin up spikes for entire duration
    spike_rate_ = histc(datarun.spikes{cell_index},datarun.triggers(1):refresh_time:datarun.stimulus.java_movie.size*refresh_time);

    % take spikes just from the relevant subset of stimulus
    spike_rate_ = spike_rate_(start_stim:end_stim);

    % translate to spike times (with duplicates for multiple spikes per time bin)
    spike_times_ = [];
    for nn = 1:max(spike_rate_)
        spike_times_ = [spike_times_; find( spike_rate_ > (nn-1) )];
    end

    
    % put into storage variables
    spike_rate(cc,:) = spike_rate_;
    spike_times{cc} = sort(spike_times_);
    
end

%%
tic
fprintf('Computing cone input in frames %d to %d...',start_stim,end_stim)
T=text_waitbar;
start_time_ = clock; % note when it started

% initialize storage variables
cone_inputs = zeros(length(stims),size(Wc,2));


% cycle through each stimulus
for ss = 1:length(stims)

    T=text_waitbar(T,ss/length(stims) - 0.01);

    % note which stim
    this_stim = stims(ss);

    % get new frame
    STAFrame = datarun.stimulus.java_movie.getFrame(this_stim-1);
    new_frame = permute(reshape(STAFrame.getBuffer,3,field_width,field_height),[3 2 1]) - .5;

    %new_frame = new_frame - 0.5;
    new_frame = reshape(new_frame,[],1);

    % convert to cone space
    cone_inputs(ss,:) = full(Wc'*double(new_frame));

end
toc

%%
template_length = 6;

cell_id = 36; 

temp_cell_index = get_cell_indices(datarun, cell_id);

% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, cell_id,...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0);
                                        
temp_cone_indices = find(selection);

% get time course from cell of interest
temp_tc = datarun.stas.time_courses{temp_cell_index};
num_frames = size(temp_tc,1);

if size(temp_tc,2) == 3
    temp_tc = temp_tc(num_frames - template_length+1:num_frames,:);
    temp_tc = sum(temp_tc, 2);
elseif size(temp_tc,2) ==1
    % cut out region of interest from TC
    temp_tc = temp_tc(num_frames - template_length+1:num_frames);
else
    error('number of color channels does not equal either 1 or 3')
end

temp_tc = temp_tc - mean(temp_tc);
temp_tc = temp_tc ./ norm(temp_tc);

% flip the impulse response to estimate the temporal filter
reverse_indices = template_length:-1:1;
impulse_filter = temp_tc(reverse_indices);
%impulse_filter = temp_tc;

filtered_cone_inputs = filter(impulse_filter, 1, cone_inputs(:,temp_cone_indices));

   

%-----------------------------------------
% COMSTRUCT THE STE AND COMPUTE STA AND STC

% get the covariance matrix of the stimulus (will be subtracted from the STC)
stimulus_cov = cov(filtered_cone_inputs);

% CONSTRUCT THE SPIKE-TRIGGERED ENSEMBLE
STE = filtered_cone_inputs(spike_times{temp_cell_index},:);

% get image frame info
cell_com = datarun.cones.rf_fits{temp_cell_index}.center;
window = 30;
frame_begin_x = cell_com(2) - window;
frame_end_x = cell_com(2) + window;
frame_begin_y = cell_com(1) - window;
frame_end_y = cell_com(1) + window;

% --- COMPUTE STA ---
temp_polarity = datarun.stas.polarities{temp_cell_index};
STA = mean(STE);
STA = STA ./ max(STA) ./2;
STA = STA * temp_polarity;

matrix_scale_factor = 1;

% --- PLOT STA ---
figure(101)
image_STA = Wc(:,temp_cone_indices) * STA';
imagesc(norm_image(matrix_scaled_up(reshape(image_STA, [320 320 3]),matrix_scale_factor)))
axis(matrix_scale_factor*[frame_begin_y frame_end_y frame_begin_x frame_end_x])
axis off
axis square
%print(101, '~/Desktop/sta.pdf', '-dpdf')

% --- COMPUTE STC ---
STA_ = mean(STE)';
z_STE = STE';
z_STE = z_STE - STA_*(STA_'*z_STE)./ norm(STA_)^2;  % project out the mean
%STE = STE - (STE*STA_') * STA_ ./ norm(STA_)^2;  % project out the mean
z_filtered_cone_inputs = filtered_cone_inputs';
z_filtered_cone_inputs = z_filtered_cone_inputs - STA_*(STA_'*z_filtered_cone_inputs)./ norm(STA_)^2;  % project out the mean

STC = cov(z_STE') - cov(z_filtered_cone_inputs');

% --- FACTORIZE THE STC ---
[PCs, eig_vals_matrix] = eig(STC);
eig_vals = diag(eig_vals_matrix);
[eig_vals, sorted_eig_indices] = sort(eig_vals, 'descend');
PCs = PCs(:, sorted_eig_indices);

% get indices to meaningful dimensions
keep_indices = find(abs(eig_vals) > 1e-10); % note this number might need to be changed depending on data
% store eigvals;
eig_vals = eig_vals(keep_indices); % ./ abs(sum(s_eig_vals));
PCs = PCs(:,keep_indices);

% --- PLOT EIGENVECTORS OF THE STC with increased variance---
num_dims = 6;
for dm = 1:num_dims
    image_PC = Wc(:,temp_cone_indices) * PCs(:,dm);
    figure(dm)
    image(norm_image(matrix_scaled_up(reshape(image_PC, [320,320,3]),matrix_scale_factor)))
    hold on
    axis(matrix_scale_factor*[frame_begin_y frame_end_y frame_begin_x frame_end_x])
    axis off
    axis square
    %save_path = ['~/Desktop/PC',num2str(dm),'.pdf'];
    %print(dm, save_path, '-dpdf')

end
% --- PLOT THE SORTED EIGENVALUE SPECTRUM ---
figure(num_dims + 1); clf;
plot(eig_vals, 'ko')


[junk,xi,yi]=roipoly;


subunit_cone_indices = inpolygon(datarun.cones.centers(:,1),datarun.cones.centers(:,2),xi,yi);
subunit_cones_centers = datarun.cones.centers(subunit_cone_indices,:);
cone_labels = 1:1:length(datarun.cones.types);
subunit_cone_labels = cone_labels(subunit_cone_indices)

% plot subunit cones
plot(datarun.cones.centers(subunit_cone_indices,1), datarun.cones.centers(subunit_cone_indices,2), 'ro')


subunit_linear_weights = datarun.cones.weights(subunit_cone_indices, temp_cell_index);

% find the complimentary set of cones with closely matched weights
substitute_cone_indices = zeros(1,length(subunit_linear_weights));
substitute_weights = substitute_cone_indices;
for cn = 1:length(subunit_linear_weights);
    flag = 2;
    diff_weights = abs(datarun.cones.weights(temp_cone_indices,temp_cell_index) - subunit_linear_weights(cn));
    [sorted_weights, sorted_indices] = sort(diff_weights, 'ascend');
    
    while ismember(sorted_indices(flag), substitute_cone_indices) 
        flag = flag + 1;
    end
    
    substitute_cone_indices(cn) = temp_cone_indices(sorted_indices(flag));
    substitute_weights(cn) = datarun.cones.weights(temp_cone_indices(sorted_indices(flag)),temp_cell_index);
end

substitute_cone_centers = datarun.cones.centers(substitute_cone_indices,:);
% plot substitute_cones
plot(substitute_cone_centers(:,1), substitute_cone_centers(:,2), 'g.')


%% make cone maps

centers = datarun.cones.centers;
% an alternative tesselation
clear subunit_cone_map scattered_cone_map
[V,C] = voronoin(centers);
width = datarun.stimulus.field_width;
height = datarun.stimulus.field_height;



subunit_cone_map = zeros(width, height);
verbose = true;
figure(1); clf; axis ij
new_counter = 0;
for cn = 1:length(C)
    if ~isempty(find(subunit_cone_labels == cn,1))
        clear temp_mask
        if all(C{cn} ~=1)
            temp_mask = poly2mask(V(C{cn},1),V(C{cn},2), height, width);
            if length(find(temp_mask ==1)) <= 40
                if verbose
                    patch(V(C{cn},1),V(C{cn},2), 'r')
                    axis([1 320 1 320])
                    axis square
                    drawnow
                    pause(0.01)
                end
                subunit_cone_map = subunit_cone_map + temp_mask;
                if new_counter == 2
                    continue
                end
            end
        end
    end
end
subunit_cone_map = subunit_cone_map';

    
scattered_cone_map = zeros(width, height);
verbose = true;
figure(1); clf; axis ij
new_counter = 0;
for cn = 1:length(C)
    if ~isempty(find(substitute_cone_indices == cn,1))
        clear temp_mask
        if all(C{cn} ~=1)
            temp_mask = poly2mask(V(C{cn},1),V(C{cn},2), height, width);
            if length(find(temp_mask ==1)) <= 40
                if verbose
                    patch(V(C{cn},1),V(C{cn},2), 'r')
                    axis([1 320 1 320])
                    axis square
                    drawnow
                    pause(0.01)
                end
                scattered_cone_map = scattered_cone_map + temp_mask;
                if new_counter == 2
                    continue
                end
            end
        end
    end
end
scattered_cone_map = scattered_cone_map';
    
%% Careful - creating concatinated maps
if ~exist('all_subunit_map', 'var')
    all_subunit_map = zeros(width, height)';
    all_scattered_map = zeros(width, height)';
end

all_subunit_map = all_subunit_map + subunit_cone_map;
all_scattered_map = all_scattered_map + scattered_cone_map;
    
figure(51)   
imagesc(all_subunit_map')
colormap gray
figure(52)
imagesc(all_scattered_map')
colormap gray

write_path = '/snle/acquisition/maps/single_cone/subunits';
dlmwrite([write_path,'map-0000', '.txt'], all_subunit_map, 'delimiter', '\t', 'newline', 'pc')
dlmwrite([write_path,'map-0001', '.txt'], all_scattered_map, 'delimiter', '\t', 'newline', 'pc')


    
    
    