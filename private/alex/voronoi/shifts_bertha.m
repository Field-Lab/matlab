datarun = load_data('/Volumes/Analysis/2015-10-06-2/d00-13-norefit/data008/data008');
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);
[inputs, refresh, duration] = get_wn_movie_ath(datarun, 'BW-2-8-0.48-11111-160x160.xml');

parameters.seed = 11111;
stimulus.stixel_width = 2;
stimulus.stixel_height = 2;

field_width = 160;
field_height = 160;

stimulus.rng_init.state = parameters.seed;
stimulus.rng_init.seed = parameters.seed;
stimulus.rng_init.state = Init_RNG_JavaStyle(stimulus.rng_init.seed);
stimulus.jitter.state = stimulus.rng_init.state;

shifts = zeros(2,duration);
tic
for i=1:duration
    if mod(i, 100)==0
        i
    end
    for j=1:field_width*field_height;
        random_uint16(stimulus.jitter.state);
    end
    jitterX = mod(double(random_uint16(stimulus.jitter.state)), stimulus.stixel_width) - stimulus.stixel_width/2;
    jitterY = mod(double(random_uint16(stimulus.jitter.state)), stimulus.stixel_height) - stimulus.stixel_height/2;
    shifts(1,i) = jitterX;
    shifts(2,i) = jitterY;
end
toc

save('/Volumes/Analysis/2015-10-06-2/jitter/shifts.mat', 'shifts')