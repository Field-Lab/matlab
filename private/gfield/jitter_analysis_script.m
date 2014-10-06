
% parameters for cone size and field size
cone_parameters.x_size = 40;
cone_parameters.y_size = 40;
cone_parameters.center_radius = 0.6;
cone_parameters.effective_radius = 3;
cone_parameters.sparse = false;

% list of cone locations
location_list = [35.2 42.6; 38.3 39.2; 40.1 36.7; 36.5 36; 39.6 32.8]-20;

% spatial profile of the RGC RF
rgc_rf = zeros(cone_parameters.x_size, cone_parameters.y_size);
for cn = 1:size(location_list,1)
    input_params = cone_parameters;
    input_params.center = location_list(cn,:);
    cone_rf = make_gaussian(input_params);
    rgc_rf = rgc_rf + cone_rf;
end
%rgc_rf = (rgc_rf *0.5) + (ones(cone_parameters.x_size, cone_parameters.y_size) * 0.5);

% plot the RF of the RGC
image(norm_image(rgc_rf))


% make time course
time_course = [0 0.2 0.6 1 0.4 0 -0.2 -0.15 -0.1 -0.05 -0.01 0 0 0]; 
plot(time_course)

movie_length = 1000;

white_movie = zeros(cone_parameters.x_size,cone_parameters.y_size,movie_length); 
subtract_factor = ones(cone_parameters.x_size, cone_parameters.y_size);
for fm = 1:movie_length;
    tmp_frame = (randi([0 1], 5, 5) * 2) - 1;
    white_movie(:,:,fm) = matrix_scaled_up(tmp_frame, 8);
end


% play movie
for fm = 1:100
    tmp= squeeze(white_movie(:,:,fm));
    imagesc(tmp)
    pause(0.3)
end


%% Compute spike train
gen_sig = zeros(1,movie_length);
for fm = 1:movie_length
    gen_sig(fm) = dot(reshape(rgc_rf, 1, []), reshape(squeeze(white_movie(:,:,fm)), 1, []));
end

gen_sig = filter(time_course, 1, gen_sig);

% normalize the generator distribution by the std
gen_sig = gen_sig ./ std(gen_sig);

NL_output = exp(gen_sig - 2);

figure(2); clf;
plot(gen_sig, 'r')
hold on
plot(NL_output, 'k')


spike_train = poissrnd(NL_output);
plot(spike_train)

%% Compute STA


