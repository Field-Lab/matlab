parameters.seed = 11111;
stimulus.stixel_width = 5;
stimulus.stixel_height = 5;

stimulus.rng_init.state = parameters.seed;
stimulus.rng_init.seed = parameters.seed;
stimulus.rng_init.state = Init_RNG_JavaStyle(stimulus.rng_init.seed);
stimulus.jitter.state = stimulus.rng_init.state;

shifts = zeros(2,72067);

for i=1:72067
    jitterX = mod(double(random_uint16(stimulus.jitter.state)), stimulus.stixel_width) - floor(stimulus.stixel_width/2);
    jitterY = mod(double(random_uint16(stimulus.jitter.state)), stimulus.stixel_height) - floor(stimulus.stixel_height/2);
    shifts(1,i) = jitterX;
    shifts(2,i) = jitterY;
end

save('/Volumes/Analysis/2008-04-30-2/jitter/shifts.mat', 'shifts')