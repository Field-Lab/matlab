% kiwi
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-05-13-3/data006/data006';
path_and_name{1,2} = 'kiwi';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% blueberry
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';
path_and_name{1,2} = 'blueberry';

% peach
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-0s-2400s/data001/data001';
path_and_name{1,2} = 'peach';

% grapes
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2007-03-27-2/data014/data014/data014';
path_and_name{1,2} = 'grapes';

% apricot
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2009-04-13-5/data005/data005';
path_and_name{1,2} = 'apricot';

% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = 'plantain';

path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2007-03-27-2/data014/data014';
path_and_name{1,2} = 'grapes'

% butterfly
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-12-12-1/data005/data005/data005';
path_and_name{1,2} = 'butterfly'

% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);

% get sta summaries 
datarun = get_sta_summaries(datarun, {3,4}, 'keep_stas', false);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_rf_portraits(datarun, {3}, 'plot_radius', 18)
plot_rf_portraits(datarun, {4}, 'plot_radius', 18)

% on cells
% from kiwi: 
on_midget_one = 1326;
on_midget_two = 7383;
plot_rf_portraits(datarun, on_midget_one, 'plot_radius', 18, 'scale_factor', 8, 'cones', true, 'cone_size', 24)
print(1, '~/Desktop/on-midget-one','-dpdf')
plot_rf_portraits(datarun, on_midget_two, 'plot_radius', 18, 'scale_factor', 8)
print(2, '~/Desktop/on-midget-two','-dpdf')

% off cells
% from kiwi: 7009
off_midget_one = 4222;
off_midget_two = 7009;
plot_rf_portraits(datarun, off_midget_one, 'plot_radius', 18, 'scale_factor', 8, 'cones', true)
print(3, '~/Desktop/off-midget-one','-dpdf')
plot_rf_portraits(datarun, off_midget_two, 'plot_radius', 18, 'scale_factor', 8, 'cones', true, 'cone_size', 24)
print(4, '~/Desktop/off-midget-two','-dpdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ON midget TCs for the L and M cone of example cell one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_rf(datarun, on_midget_one)

m_one = [145, 148];
m_two = [145, 142];
l_one = [149, 143];
l_two = [151, 139];


plot_pie_from_cone(datarun, on_midget_one, m_one, 'fig_or_axes', 1)
print(1, '~/Desktop/m-cone-one.pdf', '-dpdf')

plot_pie_from_cone(datarun, on_midget_one, m_two, 'fig_or_axes', 2)
print(2, '~/Desktop/m-cone-two.pdf', '-dpdf')

plot_pie_from_cone(datarun, on_midget_one, l_one, 'fig_or_axes', 3)
print(3, '~/Desktop/l-cone-one.pdf', '-dpdf')

plot_pie_from_cone(datarun, on_midget_one, l_two, 'fig_or_axes', 4)
print(4, '~/Desktop/l-cone-two.pdf', '-dpdf')


% for bars and time courses (not obsolete!!!!)
% plot time courses of stixels of interest
%on_mid_one_index = get_cell_indices(datarun, on_midget_one);
%datarun = load_sta(datarun, 'load_sta', on_midget_one);
%sta_size = size(datarun.stas.stas{on_mid_one_index});

% figure(1)
% temp_sig_stix = false(sta_size(1:2));
% temp_sig_stix(145,148) = true;
% temp_tc = time_course_from_sta(datarun.stas.stas{on_mid_one_index}, temp_sig_stix);
% max_green = max(temp_tc(:,2));
% clf
% hold on
% plot(temp_tc(:,1) ./ max_green, 'r')
% plot(temp_tc(:,2) ./ max_green, 'g')
% plot(temp_tc(:,3) ./ max_green, 'b')
% axis([1 6 -0.5 1])
% hold off
% print(1, '~/Desktop/m-cone.pdf', '-dpdf')
% 
% figure(2)
% temp_sig_stix = false(sta_size(1:2));
% temp_sig_stix(149,143) = true;
% temp_tc = time_course_from_sta(datarun.stas.stas{on_mid_one_index}, temp_sig_stix);
% max_green = max(temp_tc(:,2));
% clf
% hold on
% plot(temp_tc(:,1) ./ max_green, 'r')
% plot(temp_tc(:,2) ./ max_green, 'g')
% plot(temp_tc(:,3) ./ max_green, 'b')
% axis([1 6 -0.5 1])
% hold off
% print(2, '~/Desktop/1-cone.pdf', '-dpdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Off midget TCs for the L and M cone of example cell one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_rf(datarun, off_midget_one, 'polarity', true)

% assign stixels of interest
% m cones
m_one = [129, 185];
m_two = [131, 182];
l_one = [131, 189];
l_two = [133, 186];


plot_pie_from_cone(datarun, off_midget_one, m_one, 'fig_or_axes', 1)
print(1, '~/Desktop/m-cone-one.pdf', '-dpdf')

plot_pie_from_cone(datarun, off_midget_one, m_two, 'fig_or_axes', 2)
print(2, '~/Desktop/m-cone-two.pdf', '-dpdf')

plot_pie_from_cone(datarun, off_midget_one, l_one, 'fig_or_axes', 3)
print(3, '~/Desktop/l-cone-one.pdf', '-dpdf')

plot_pie_from_cone(datarun, off_midget_one, l_two, 'fig_or_axes', 4)
print(4, '~/Desktop/l-cone-two.pdf', '-dpdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Classification plots

datarun = load_data('2008-05-13-3','rf-6-kiwi');
%datarun = load_data('2009-04-13-5', 'rf-5-apricot');
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = get_sta_summaries(datarun, {1,2,3,4,5}, 'keep_stas', false);

% load cone weights matrix
load([single_cone_path datarun.names.nickname '/Wc.mat'])
connections_path = [single_cone_path datarun.names.nickname '/connections.txt'];
connections = dlmread(connections_path);

[num_cones, num_rgcs] = size(connections);

trusted_rgc_types = {1,2,3,4,5};
rgc_indices = get_cell_indices(datarun, trusted_rgc_types);
trusted_connections = connections;

significance_threshold = 4.0;  % units of sigmas
false_cone_counter = 0;
kept_cone_counter = 0;
clear rejected_cones kept_cone_indices norm_rgbs

for cn = 1:num_cones
    % cone vector
    cn_vector = reshape(full(Wc(:,cn)),[],3);
    % force each 'color' vector to have unit length.
    for clm = 1:3
        cn_vector(:,clm) = cn_vector(:,clm) ./ norm(cn_vector(:,clm));
    end
    
    % compute threshold for a RGC to contribute to cone typing from the 
    %distribution of connections between cone and all rgcs.
    sampling_thresh = std(trusted_connections(cn,:)) * significance_threshold;
    
    % get indices to these rgcs with strongest connections to cone
    rgcs_of_interest_center = find(trusted_connections(cn,:) > sampling_thresh);
    rgcs_of_interest_surround = find(trusted_connections(cn,:) < -sampling_thresh);
    rgcs_of_interest = union(rgcs_of_interest_center, rgcs_of_interest_surround);

    %%% troubleshoot
    tmp_rgcs = datarun.cell_ids(rgc_indices(rgcs_of_interest));
    

    if isempty(rgcs_of_interest)
        false_cone_counter = false_cone_counter + 1
        rejected_cones(false_cone_counter) = cn
    else
        % intialize some variables that will be used during loop
        kept_cone_counter = kept_cone_counter + 1;
        kept_cone_indices(kept_cone_counter) = cn;
        rgb_estimates = zeros(3,length(rgcs_of_interest));
        variances = zeros(1,length(rgcs_of_interest));

        % for each significant rgc, get rgb triplet associated with cone
        % add together rgb estimates in units of SD (~optimal)
        for rgc = 1:length(rgcs_of_interest)
            temp_weights = trusted_connections(:,rgcs_of_interest(rgc)); % weights rgc
            temp_weights = temp_weights ./ robust_std(temp_weights); % normalize weights by robust SD

            variances(rgc) = 1/temp_weights(cn);  

            temp_rf = reshape(datarun.stas.rfs{rgc_indices(rgcs_of_interest(rgc))}, [],3);
            temp_rgb = diag(temp_rf' * cn_vector)';
            rgb_estimates(:,rgc) = temp_rgb;        
        end
        rgb = variance_weighted_average(rgb_estimates, variances);
        if isnan(rgb)
            cn
        end

        % store rgb values in datarun and make sure they are on the unit sphere
        norm_rgbs(kept_cone_counter,:) = rgb ./ norm(rgb);
    end
end
neg_indices = find(norm_rgbs(:,2) < 0);
norm_rgbs(neg_indices,:) = norm_rgbs(neg_indices,:) * -1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datarun = import_single_cone_data(datarun);
datarun.piece.rig = 'A';
datarun.piece.optical_path_direction = 'below';

%plot_cone_classification(datarun, 'foa_3d', 1, 'foa_pie', 2, 'foa_2d', 3, 'foa_hist', 4)
%print(1, '~/Desktop/3D-scatter.pdf', '-dpdf')
%print(4, '~/Desktop/classification_hist.pdf', '-dpdf')

%%% alternative classification
s_cones = find(datarun.cones.types(kept_cone_indices) == 'S');
l_cones = find(datarun.cones.types(kept_cone_indices) == 'L');
m_cones = find(datarun.cones.types(kept_cone_indices) == 'M');
u_cones = find(datarun.cones.types(kept_cone_indices) == 'U');
y_cones = find(datarun.cones.types(kept_cone_indices) ~= 'S');

rgb_expected = cone_rgb_expected(datarun);

% flip fucked S cones
neg_indices = find(norm_rgbs(s_cones,3) < 0);
norm_rgbs(s_cones(neg_indices),:) = norm_rgbs(s_cones(neg_indices),:) * -1;

figure(1)
clf

[X, Y, Z] = sphere(52);
% removing parts of the sphere
X_pos = find(X < 0);
Y_pos = find(Y < 0);
Z_pos = find(Z < 0);
%X(X_pos) = NaN;
Y(Y_pos) = NaN;
Z(Z_pos) = NaN;

a_data = zeros(53, 53);
a_data(:,:) = 0.3;
mesh(X, Y, Z, 'AlphaDataMapping', 'none','FaceAlpha', 0.3)
%surf(X, Y, Z, 'AlphaDataMapping', 'none','FaceAlpha', 0.3)
%shading interp
null_colormap = [0 0 0; 0 0 0; 0 0 0];
hidden off
colormap(null_colormap)
hold on
axis square

plot3(norm_rgbs(s_cones,1), norm_rgbs(s_cones,2), norm_rgbs(s_cones,3), '.', 'Color', [1 1 1])
plot3(norm_rgbs(l_cones,1), norm_rgbs(l_cones,2), norm_rgbs(l_cones,3), '.', 'Color', [1 1 1])
plot3(norm_rgbs(m_cones,1), norm_rgbs(m_cones,2), norm_rgbs(m_cones,3), '.', 'Color', [1 1 1])

plot3(norm_rgbs(s_cones,1), norm_rgbs(s_cones,2), norm_rgbs(s_cones,3), '.', 'Color', [0.4 0.4 1])
plot3(norm_rgbs(l_cones,1), norm_rgbs(l_cones,2), norm_rgbs(l_cones,3), '.', 'Color', [1 0.4 0.4])
plot3(norm_rgbs(m_cones,1), norm_rgbs(m_cones,2), norm_rgbs(m_cones,3), '.', 'Color', [0.4 1 0.4])

%plot3(norm_rgbs(m_cones(sorted_indices(nth_furthest_cone)),1), norm_rgbs(m_cones(sorted_indices(nth_furthest_cone)),2), norm_rgbs(m_cones(sorted_indices(nth_furthest_cone)),3), 'k.')

%plot gun axes
plot3([-2 2], [0 0], [0 0], 'k')
plot3([0 0], [-2 2], [0 0], 'k')
plot3([0 0], [0 0], [-2 2], 'k')

% plot cone axes
plot3([0 2*rgb_expected.L(1)], [0 2*rgb_expected.L(2)], [0 2*rgb_expected.L(3)], 'Color', [0.75 0 0], 'LineWidth', 2)
plot3([0 2*rgb_expected.M(1)], [0 2*rgb_expected.M(2)], [0 2*rgb_expected.M(3)], 'Color', [0 0.75 0], 'LineWidth', 2)
plot3([0 2*rgb_expected.S(1)], [0 2*rgb_expected.S(2)], [0 2*rgb_expected.S(3)], 'Color', [0 0 0.75], 'LineWidth', 2)

view([-207, 30])

hold off
print(1, '~/Desktop/classification.pdf', '-dpdf')

%%%%%%%%%%%%%%%%%%%%%
% make L/M histogram

mean_L = mean(norm_rgbs(l_cones,:))
mean_M = mean(norm_rgbs(m_cones,:))
mean_S = robust_mean(norm_rgbs(s_cones,:))

discrim_axis = mean_L - mean_M;
proj_L = norm_rgbs(l_cones,:) * discrim_axis';
proj_M = norm_rgbs(m_cones,:) * discrim_axis';

% make histogram 
hist_bins = [-0.15:0.005:0.15];
[l_hist, hist_bins] = hist(proj_L, hist_bins);
[m_hist, hist_bins] = hist(proj_M, hist_bins);

figure(9); clf;
bar(hist_bins,l_hist, 'r')
hold on
bar(hist_bins, m_hist, 'g')
hold off
axis([-0.125 0.125 0 160])
print(9, '~/Desktop/lm_histogram.pdf', '-dpdf')


num_L = length(l_cones);
num_M = length(m_cones);
groups = zeros(num_L+num_M,1);
groups(1:num_L) = 1;
[classed, est_error] = classify(rand(50,1), [proj_L; proj_M],groups);


%%%%%%%%%%%%%%%%%%%%%
% troubleshoot

m_mean = mean(norm_rgbs(m_cones,:));
m_dists = ipdm(norm_rgbs(m_cones,:), m_mean);

[m_dists_sorted, sorted_indices] = sort(m_dists, 'descend');

kept_cone_indices(m_cones(sorted_indices(1)))


%%%%%%%%%%%%%%%%%%%%%


% make 2D scatter plot
figure(1)
clf
loglog(rgbs(s_cones,3)./rgbs(s_cones,2), rgbs(s_cones,1)./rgbs(s_cones,2), 'b.','MarkerSize', 8)
hold on
loglog(rgbs(l_cones,3)./rgbs(l_cones,2), rgbs(l_cones,1)./rgbs(l_cones,2), 'r.', 'MarkerSize', 8)
loglog(rgbs(m_cones,3)./rgbs(m_cones,2), rgbs(m_cones,1)./rgbs(m_cones,2), 'g.', 'MarkerSize', 8)
loglog(rgbs(u_cones,3)./rgbs(u_cones,2), rgbs(u_cones,1)./rgbs(u_cones,2), '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 8)
loglog(rgb_expected.S(3)./rgb_expected.S(2), rgb_expected.S(1)./rgb_expected.S(2), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
loglog(rgb_expected.M(3)./rgb_expected.M(2), rgb_expected.M(1)./rgb_expected.M(2), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
loglog(rgb_expected.L(3)./rgb_expected.L(2), rgb_expected.L(1)./rgb_expected.L(2), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
hold off
print(1, '~/Desktop/2D-scatter.pdf', '-dpdf')



% separate S cones from L,M and U cones
hist_bins = [-5:0.2: 10];
[hist_all, hist_bins] = hist(log(rgbs(y_cones,3)./rgbs(y_cones,2)), hist_bins);
[hist_blue, hist_bins] = hist(log(rgbs(s_cones,3)./rgbs(s_cones,2)), hist_bins);

figure(1)
clf
hold on
bar(hist_bins, hist_all, 'k')
bar(hist_bins, hist_blue, 'b')
hold off
print(1, '~/Desktop/s-cone-hist', '-dpdf')

% separate L and M cones
discriminator = (rgbs(:,1)./sum(rgbs,2)) - (rgbs(:,3) ./ sum(rgbs,2));

[hist_test, hist_bins] = hist(discriminator, 20);
plot(hist_bins, hist_test)


[hist_green, hist_bins] = hist(log(rgbs(s_cones,3)./rgbs(s_cones,2)), hist_bins);

figure(2)
clf
hold on
bar(hist_bins, hist_all, 'k')
bar(hist_bins, hist_blue, 'b')
hold off
print(2, '~/Desktop/l-m-cone-hist.pdf', '-dpdf')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cone mosaic plot
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);
plot_cone_mosaic(datarun, 'fig_or_axes', 1, 'cone_size', 6)
print(1, '~/Desktop/cone-mosaic.pdf', '-dpdf')

plot_cell_sampling(datarun, [1326 4222], 'fig_or_axes', 1, 'plot_cones', true, 'thresh', 0.1, 'contiguity', true, 'scale', 2.5, 'cone_size', 6, 'label', true)

plot_cell_sampling(datarun, {4}, 'fig_or_axes', 2, 'plot_cones', true, 'thresh', 0.05, 'contiguity', true, 'scale', 2.5, 'cone_size', 6, 'label', true)

plot_rf_portraits(datarun, [1054 589 1067 2809 650 1321],'plot_radius', 18)
