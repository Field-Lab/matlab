function trial_dur = get_mb_trial_dur(datarun)

display_width = 800; 
display_height = 600;
refresh_rate = 60.35;
delta = datarun.stimulus.params.DELTA;
bar_width = datarun.stimulus.params.BAR_WIDTH;

for dt = 1:length(delta)
    trial_dur(dt) = (sqrt(display_width^2+display_height^2)+bar_width)/delta(dt)/refresh_rate;
end  