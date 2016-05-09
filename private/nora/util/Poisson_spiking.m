function logical_sim = Poisson_spiking(firing_rate, trials, bins_per_frame, monitor_refresh)

frames = length(firing_rate);
bins = length(firing_rate)*bins_per_frame;
logical_sim = zeros(trials, bins);
bindur = 1/(monitor_refresh*bins_per_frame);
for i_trial = 1 : trials
    bin = 1;
    binary_simulation = zeros(1,bins);
    for i = 1 : frames;
        for bins_in_frame = 1:10
            roll = rand(1);
            if roll >  exp(-bindur*firing_rate(i));
                binary_simulation(bin)= 1;
            end
            bin = bin+1;
        end
    end
    logical_sim(i_trial,:) = binary_simulation ;
end
end

