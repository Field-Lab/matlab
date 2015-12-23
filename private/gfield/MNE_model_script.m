% Get the information for a single cell
[datarun, cone_info] = load_data_and_cones('plantain', 'sta_summaries', false);

[datarun, cone_info] = load_data_and_cones('apple', 'sta_summaries', false);


datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);
datarun = set_polarities(datarun);
plot_matrix_scale_factor = 1;
Nlags = 1; % == 1, corresponds to a spatial rf
save_path = '~/Analysis/Minimal-Models/plantain/';

height = datarun.stimulus.field_height;
width = datarun.stimulus.field_width;

order   = 2;   % order of MNE model to fit
fittype = 0;   % 0 for regular fitting, 1 for random fitting
jack = 4;   % # jackknives to run (also determines the size of each jackknives)
Nlags   = 1;   % # time lags to use in spatiotemporal models (=1 is just a spatial model)



%% Select cell of interest and get connected cones 

%plantain
cell_id = 497; % on midget plantain
cell_id = 696; % weak on parasol in plantain

cell_id = 947; % off parasol from apple

cell_index = get_cell_indices(datarun, cell_id);
[mosaic_weights, selection, extras] = select_cone_weights(datarun, cell_id,...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0);
cone_ids = find(selection);                                        

% get TC for cell of interest and filter
datarun = get_sta_summaries(datarun, cell_id, 'keep_stas', false);
temp_tc = datarun.stas.time_courses{cell_index};
temp_tc = temp_tc(:,2)
temp_tc = temp_tc - mean(temp_tc);
temp_tc = temp_tc ./ norm(temp_tc);

% flip the impulse response to estimate the temporal filter
reverse_indices = length(temp_tc):-1:1;
impulse_filter = temp_tc(reverse_indices);

raw_cone_stimulus = cone_info.cone_inputs(:,cone_ids);


% filter the cone projected stimulus by the RGC time course
stimulus = filter(impulse_filter, 1, raw_cone_stimulus);


% get center of mass for this cell
window_size = 25;
rf_center = datarun.stas.rf_coms{cell_index};
bg_window_y = floor(rf_center(2) - window_size);
ed_window_y = floor(rf_center(2) + window_size);
bg_window_x = floor(rf_center(1) - window_size);
ed_window_x = floor(rf_center(1) + window_size);

%% fit with maximum entroy model
[Nsamples, Ndim] = size(stimulus);

% The stimulus units is changed to be z-score
% (zero mean, normalized by standard deviation),
% but could be put into your own unit system (or 
% maybe your stimulus file is already saved
% in the proper units).
stimulus = stimulus-repmat(mean(stimulus),[Nsamples,1]);
stimulus = stimulus./repmat(std(stimulus),[Nsamples,1]);


njacks = 4;
test_fit = zeros(njacks, 1+length(cone_ids)+(length(cone_ids).^2));
for jack = 1:njacks
    Nsamples = Nsamples/4;
    resp = cone_info.spike_rate(cell_index,:)';
    stim = stimulus;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% break into training and test sets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ntest = floor(Nsamples/4);  % could be changed to do different # of jackknives
    teststim = stim(1+(jack-1)*Ntest:jack*Ntest,:);
    testresp = resp(1+(jack-1)*Ntest:jack*Ntest);
    stim(1+(jack-1)*Ntest:jack*Ntest,:) = [];
    resp(1+(jack-1)*Ntest:jack*Ntest) = [];


    disp('Starting optimization');
    tic
    celltype = '';  % ignore this
    temp_fit = MNEfit(stim, resp, teststim, testresp, celltype, cell_id, jack, order, save_path, Ndim, fittype);
    test_fit(jack,:) = temp_fit;
    
    disp(['Optimization took ' num2str(toc/60) ' minutes']);
end

%% load results from fitting max entropy model 
njacks = 4;
file_path = '~/Analysis/Minimal-Models/plantain/';
cell_id = 947;

for jk = 1:njacks
    load_path = [file_path,'ModelCell_',num2str(cell_id),'_MNE2_jack_',num2str(jk),'.mat'];
    temp_fit = load(load_path);
    test_fit(jk,:) = temp_fit.pfinal;
end



%%

a = mean(test_fit(1:njacks,1));
h = mean(test_fit(1:njacks,2:length(cone_ids)+1));
J = mean(test_fit(1:njacks,length(cone_ids)+2:end));

J = reshape(J, length(cone_ids),length(cone_ids));

plot_matrix_scale_factor = 5;

% --- Prepare figure ---
figure(10); clf;
panel_num= 12;
px = ceil(sqrt(panel_num));
py = ceil(panel_num)/px;
plot_axes = subplot_axes(10,[0 0 1 .95],0.05,0.15,px,py);



% --- PLOT linear term h ---
cone_STA = cone_info.Wc(:,cone_ids) * -h';
cone_rf = reshape(cone_STA, [height, width, 3]);
cone_rf = cone_rf(bg_window_y:ed_window_y,bg_window_x:ed_window_x,:);
axes(plot_axes{1})
imagesc(norm_image(matrix_scaled_up(cone_rf, plot_matrix_scale_factor)))
axis off; axis image
title('linear term')


% get the eigenvals and vectors for the quadratic term
[V, D] = eig(J);
[sorted_eigvals, sorted_indices] = sort(diag(D), 'ascend');
sorted_eigenvecs = V(:,sorted_indices);

axes(plot_axes{2})
plot(sorted_eigvals, 'k.')
title('eigenvalue spectrum')

% --- PLOT EIGENVECTORS OF THE STC with increased variance---
num_dims = 10
for dm = 1:num_dims
    %figure(20+dm)
    axes(plot_axes{dm + 2})
    image_PC = cone_info.Wc(:,cone_ids) * sorted_eigenvecs(:,dm);
    image_PC = reshape(image_PC, [height, width, 3]);
    image_PC = image_PC(bg_window_y:ed_window_y,bg_window_x:ed_window_x,:);
    image(norm_image(matrix_scaled_up(image_PC, plot_matrix_scale_factor)))
    axis off; axis square
    title(['D',num2str(dm)])
    drawnow
end

print(10, '~/desktop/kernels.pdf', '-dpdf')


%% check the individual jackknives

for jk = 1:njacks
    
    a = test_fit(jk,1);
    h = test_fit(jk,2:length(cone_ids)+1);
    J = test_fit(jk,length(cone_ids)+2:end);

    J = reshape(J, length(cone_ids),length(cone_ids));

    plot_matrix_scale_factor = 5;

    % --- Prepare figure ---
    figure(jk); clf;
    px = ceil(sqrt(9));
    py = ceil(9)/px;
    plot_axes = subplot_axes(jk,[0 0 1 .95],0.05,0.15,px,py);



    % --- PLOT linear term h ---
    cone_STA = cone_info.Wc(:,cone_ids) * -h';
    cone_rf = reshape(cone_STA, [height, width, 3]);
    cone_rf = cone_rf(bg_window_y:ed_window_y,bg_window_x:ed_window_x,:);
    axes(plot_axes{1})
    imagesc(norm_image(matrix_scaled_up(cone_rf, plot_matrix_scale_factor)))
    axis off; axis image
    title('linear term')


    % get the eigenvals and vectors for the quadratic term
    [V, D] = eig(J);
    [sorted_eigvals, sorted_indices] = sort(diag(D), 'ascend');
    sorted_eigenvecs = V(:,sorted_indices);

    % --- PLOT EIGENVECTORS OF THE STC with increased variance---
    num_dims = 25;
    for dm = 10:num_dims+10
        %figure(20+dm)
        axes(plot_axes{dm + 1})
        image_PC = cone_info.Wc(:,cone_ids) * sorted_eigenvecs(:,dm);
        image_PC = reshape(image_PC, [height, width, 3]);
        image_PC = image_PC(bg_window_y:ed_window_y,bg_window_x:ed_window_x,:);
        image(norm_image(matrix_scaled_up(image_PC, plot_matrix_scale_factor)))
        axis off; axis square
        title(['D',num2str(dm)])
        drawnow
    end

end

%%

J_prime = J' - h'*(h*J')./ norm(h)^2;  % project out the mean



Y = pdist(J,'cityblock'); 
Z = linkage(Y,'average'); 
[H, T] = dendrogram(Z);








