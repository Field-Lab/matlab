% apricot classification figure
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2009-04-13-5/data008/data008';


% load data
datarun = load_data(path_and_name{1,1});
datarun = load_neurons(datarun);
datarun = load_index(datarun);


datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);

obvius_fit_path = '/snle/lab/Experiments/Array/Analysis/2009-04-13-5/rf-8-gf/stas-008/';
datarun.names.obvius_fit_path = obvius_fit_path;

datarun = load_obvius_sta_fits(datarun);
datarun = get_sta_fits_from_obvius(datarun, {1,2,3,4,5});
datarun = load_monitor_alignment(datarun);

% get sta summaries 
datarun = get_sta_summaries(datarun, {1,2,3,4,5}, 'keep_stas', false);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% on parasol rfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% on parasol mosaic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot_rf_fit(datarun, {1}, 'fits_to_use', 'obvius','fig_or_axes', 1)
plot_rf_summaries(datarun,{1},'label',false,'array',true,'array_label',false)
axis([0 80 0 40])
print(1, '~/Desktop/on-parasol-mosaic','-deps')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% off parasol rfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% off parasol mosaic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_rf_summaries(datarun,{2},'label',false,'array',true,'array_label',false)
axis([0 80 0 40])
print(1, '~/Desktop/off-parasol-mosaic','-deps')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% on midget rfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% on midget mosaic 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_rf_summaries(datarun,{3},'label',false,'array',true,'array_label',false)
axis([0 80 0 40])
print(1, '~/Desktop/on-midget-mosaic','-deps')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% off midget rfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% off midget mosaic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_rf_summaries(datarun,{4},'label',false,'array',true,'array_label',false)
axis([0 80 0 40])
print(1, '~/Desktop/off-midget-mosaic','-deps')
print(1, '~/Desktop/off-midget-mosaic','-dpdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBC rfs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% portraits
plot_rf_summaries(datarun,{5},'label',false,'array',true,'array_label',false)
axis([0 80 0 40])
print(1, '~/Desktop/sbc-mosaic','-deps')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Classification Fig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCA on time courses

cell_types = {1,2,3,4,5};
temp_cell_indices = get_cell_indices(datarun, cell_types);

sbc_indices = get_cell_indices(datarun, {5});
sbc_pointers = ismember(temp_cell_indices,sbc_indices);
on_m_indices = get_cell_indices(datarun, {3});
on_m_pointers = ismember(temp_cell_indices, on_m_indices);
off_m_indices = get_cell_indices(datarun, {4});
off_m_pointers = ismember(temp_cell_indices, off_m_indices);
off_p_indices = get_cell_indices(datarun, {2});
off_p_pointers = ismember(temp_cell_indices, off_p_indices);
on_p_indices = get_cell_indices(datarun, {1});
on_p_pointers = ismember(temp_cell_indices, on_p_indices);


num_rgcs = length(temp_cell_indices);
[tc_length, num_colors] = size(datarun.stas.time_courses{temp_cell_indices(1)});

tc_ensemble = zeros(num_rgcs, (tc_length * num_colors));
radii = zeros(num_rgcs,1);
green_blue = zeros(num_rgcs,1);

for cc = 1:num_rgcs
    temp_tc = datarun.stas.time_courses{temp_cell_indices(cc)}(25:30,:);
    temp_blue_signal = sum(abs(temp_tc(:,3))) * sign(ext(temp_tc(:,3)));
    temp_green_signal = sum(abs(temp_tc(:,2))) * sign(ext(temp_tc(:,2)));
    temp_red_signal = sum(abs(temp_tc(:,1))) * sign(ext(temp_tc(:,1)));
    %green_blue(cc) = temp_blue_signal - temp_green_signal - temp_red_signal ./ (abs(temp_green_signal) + abs(temp_blue_signal) + abs(temp_red_signal));
    green_blue(cc) = temp_green_signal - temp_blue_signal;
    tc_ensemble(cc,:) = reshape(datarun.stas.time_courses{temp_cell_indices(cc)},1,[]);
    radii(cc) = get_rf_fit_radius(datarun, datarun.cell_ids(temp_cell_indices(cc)), 'fits_to_use', 'obvius', 'units', 'microns');
end
    
[PCs, weights, eigenvals] = princomp(tc_ensemble);

plot(weights(:,1), weights(:,2), 'ko');

figure(1)
plot(weights(sbc_pointers,1)*-1, radii(sbc_pointers), 'ko');
hold on
plot(weights(on_m_pointers,1)*-1, radii(on_m_pointers), 'ro');
plot(weights(off_m_pointers,1)*-1, radii(off_m_pointers), 'bo');
plot(weights(off_p_pointers,1)*-1, radii(off_p_pointers), 'co');
plot(weights(on_p_pointers,1)*-1, radii(on_p_pointers), 'mo');
hold off
axis([-0.75 0.75 0 120])
axis square
print(1, '~/Desktop/classes.pdf', '-dpdf')


plot(green_blue, weights(:,1), 'ko')
hold on
plot(green_blue(sbc_pointers), weights(sbc_pointers,1), 'ro');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time courses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the times courses for parasol cells

figure(1); hold on;
on_parasol_indices = get_cell_indices(datarun, {1});
for rgc = 1:length(on_parasol_indices)
    
    tc_matrix = datarun.stas.time_courses{on_parasol_indices(rgc)}(20:30,:);
    % normalize by the variance of the green gun
    std_norm = std(tc_matrix(:,2));
    norm_tc_matrix = tc_matrix ./ std_norm;
    
    plot(norm_tc_matrix(:,1), 'r')
    plot(norm_tc_matrix(:,2), 'g')
    plot(norm_tc_matrix(:,3), 'b')
end

figure(3); hold on;
on_midget_indices = get_cell_indices(datarun, {3});
for rgc = 1:length(on_parasol_indices)
    plot(datarun.stas.time_courses{on_midget_indices(rgc)}(:,1), 'r')
    plot(datarun.stas.time_courses{on_midget_indices(rgc)}(:,2), 'g')
    plot(datarun.stas.time_courses{on_midget_indices(rgc)}(:,3), 'b')
end












