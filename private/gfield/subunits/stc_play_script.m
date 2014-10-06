% Get the information for a single cell
[datarun, cone_info] = load_data_and_cones('apple', 'sta_summaries', false);
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);
datarun = set_polarities(datarun);


plot_matrix_scale_factor = 1;

%%%%%%%
% get index to RGC of interest
% peach
%cell_id = 3334; % a nice off parasol from peach
%cell_id = 7051;

% apple
cell_id = 947;

% plantain
cell_id = 497

%%
cell_index = get_cell_indices(datarun, cell_id);

[mosaic_weights, selection, extras] = select_cone_weights(datarun, cell_id,...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0);
                                        
cone_ids = find(selection);                                        

%%
width = datarun.stimulus.field_width;
height = datarun.stimulus.field_height;
midget_window_size = 15;
parasol_window_size = 40;
                
temp_cell_id = datarun.cell_ids(cell_index);
% -- Determine window size for plotting --
if ismember(temp_cell_id, datarun.cell_types{3}.cell_ids) || ismember(temp_cell_id, datarun.cell_types{4}.cell_ids)
    window_size = midget_window_size;
else
    window_size = parasol_window_size;
end

% --- Prepare figure ---
% figure(1); clf;
% px = ceil(sqrt(9));
% py = ceil(9)/px;
% plot_axes = subplot_axes(1,[0 0 1 .95],0.05,0.15,px,py);


% get the portrait for the RF
temp_rf = get_rf(datarun, cell_id, 'polarity', true);
temp_com = round(datarun.stas.rf_coms{cell_index,:});
bg_window_x = (temp_com(1) - window_size); %* plot_matrix_scale_factor;
ed_window_x = (temp_com(1) + window_size); %* plot_matrix_scale_factor;
bg_window_y = (temp_com(2) - window_size); %* plot_matrix_scale_factor;
ed_window_y = (temp_com(2) + window_size); %* plot_matrix_scale_factor;
% ensure that sta image can fit in window, otherwise skip cell.
if ed_window_y > height || ed_window_x > width || bg_window_x < 0 || bg_window_y < 0
    continue
end
%%

% --- GET STA AND STC IN CONE SPACE ---

% get time course from cell of interest, zero the mean and make unit
% lenth
temp_tc = cone_info.rgc_tcs(cell_index,:);
temp_tc = temp_tc - mean(temp_tc);
temp_tc = temp_tc ./ norm(temp_tc);

% flip the impulse response to estimate the temporal filter
reverse_indices = length(temp_tc):-1:1;
impulse_filter = temp_tc(reverse_indices);

% filter the cone projected stimulus by the RGC time course
filtered_cone_inputs = filter(impulse_filter, 1, cone_info.cone_inputs(:,cone_ids));



% --- COMSTRUCT THE STE AND COMPUTE STA AND STC ---


% get the covariance matrix of the stimulus (will be subtracted from the STC)
stimulus_cov = cov(filtered_cone_inputs);

% CONSTRUCT THE SPIKE-TRIGGERED ENSEMBLE
STE = filtered_cone_inputs(cone_info.spike_times{cell_index},:);

%%    



% --- COMPUTE STA ---
temp_polarity = datarun.stas.polarities{cell_index};
STA = mean(STE);
STA = STA ./ max(STA) ./2;
STA = STA * temp_polarity;


% --- PLOT STA ---
%axes(plot_axes{2})
cone_STA = cone_info.Wc(:,cone_ids) * STA';
cone_rf = reshape(cone_STA, [height, width, 3]);
cone_rf = cone_rf(bg_window_y:ed_window_y,bg_window_x:ed_window_x,:);
imagesc(norm_image(matrix_scaled_up(cone_rf, plot_matrix_scale_factor)))
axis off; axis image
title('cone space from STA')


% --- COMPUTE STC ---
STA_ = mean(STE)';
z_STE = STE';
z_STE = z_STE - STA_*(STA_'*z_STE)./ norm(STA_)^2;  % project out the mean
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

% --- PLOT THE SORTED EIGENVALUE SPECTRUM ---
plot(eig_vals, 'ko')

% --- PLOT EIGENVECTORS OF THE STC with increased variance---
num_dims = 6;
for dm = 1:num_dims
    figure(20+dm)
    image_PC = cone_info.Wc(:,cone_ids) * PCs(:,dm);
    image_PC = reshape(image_PC, [height, width, 3]);
    image_PC = image_PC(bg_window_y:ed_window_y,bg_window_x:ed_window_x,:);
    image(norm_image(matrix_scaled_up(image_PC, plot_matrix_scale_factor)))
    axis off; axis square
    title(['D',num2str(dm)])
    drawnow
end

%%

num_subunits = 7;
IDX = kmeans(PCs(:,1:6),num_subunits);

% --- PLOT EIGENVECTORS OF THE STC with increased variance---
num_dims = num_subunits;
for dm = 1:num_dims

    temp_IDs = find(IDX == dm);
    ID_indices = zeros(1,length(cone_ids));
    ID_indices(temp_IDs) = 1;
    
    subunit_image = cone_info.Wc(:,cone_ids) * ID_indices';
    figure(dm+20)
    subunit_image = reshape(subunit_image, [height, width, 3]);
    subunit_image = subunit_image(bg_window_y:ed_window_y,bg_window_x:ed_window_x,:);
    image(norm_image(subunit_image))

%   image(matrix_scaled_up(norm_image(reshape(subunit_image, [320,320,3])), matrix_scale_factor))
%   axis(matrix_scale_factor*[frame_begin_y frame_end_y frame_begin_x frame_end_x])
    axis off; axis square
%    save_path = ['~/Desktop/PC',num2str(dm),'.pdf'];
%    print(num_dims+dm, save_path, '-dpdf', '-painters')
end

%%

