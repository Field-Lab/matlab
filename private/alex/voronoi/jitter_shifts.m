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




% % % 2016-02-17-6 data026
% % % RGB-16-2-0.48-22222-119.5

parameters.seed = 22222;
stimulus.stixel_width = 16;
stimulus.stixel_height = 16;

field_width = 40;
field_height = 20;

stimulus.rng_init.state = parameters.seed;
stimulus.rng_init.seed = parameters.seed;
stimulus.rng_init.state = Init_RNG_JavaStyle(stimulus.rng_init.seed);
stimulus.jitter.state = stimulus.rng_init.state;

duration = 107551;
shifts = zeros(2,duration);
tic
for i=1:duration
    for j=1:800
        random_uint16(stimulus.jitter.state);
    end
    jitterX = mod(double(random_uint16(stimulus.jitter.state)), stimulus.stixel_width) - floor(stimulus.stixel_width/2);
    jitterY = mod(double(random_uint16(stimulus.jitter.state)), stimulus.stixel_height) - floor(stimulus.stixel_height/2);
    shifts(1,i) = jitterX;
    shifts(2,i) = jitterY;
end
toc

save('/Volumes/Analysis-1/2016-02-17-6/jitter/shifts.mat', 'shifts')

new_inputs = zeros(800,3,107301);
for j=2:3
    a = squeeze(inputs(:,j,:));
    tmp = a(:);
    adj = 0;
    for i=1:100000
        tmp(i*800+1+adj:i*800+2+adj) = -1;
        adj=adj+2;
    end
    tmp(tmp==-1) = [];
    new_inputs(:,j,:) = reshape(tmp,800,107301);
end
clear tmp a



