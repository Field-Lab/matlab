function adj_TS = find_glitches(time_stamps, triggers)
    adj_TS = time_stamps;
    TS_Start = 1;
    %block_length = frames_per_trigger/monitor_refresh;
    for i=1:(length(triggers))
        supposed_t_start=triggers(i);
        actual_t_start = time_stamps(TS_Start);
        adj_TS(TS_Start+(0:99))=adj_TS(TS_Start+(0:99))+supposed_t_start-actual_t_start;
        TS_Start = TS_Start + 100;
    end
end