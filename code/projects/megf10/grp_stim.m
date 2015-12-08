function index = grp_stim(datarun)

% index = grp_stim(datarun)
% group stimulus trials by stimulus combination, return the index
% index: M spatial period x N temporal period x R directions x S repeats
%
% xyao
% 2013-12-16

list = datarun.stimulus.trial_list;
repeat = datarun.stimulus.repetitions;
stim_comb = length(datarun.stimulus.combinations);

sp = datarun.stimulus.params.SPATIAL_PERIOD;
tp = datarun.stimulus.params.TEMPORAL_PERIOD;
dr = datarun.stimulus.params.DIRECTION;



index_temp = zeros(stim_comb, repeat);

for i = 1:stim_comb
    index_temp(i, :) = find(list == i);
end

[sp_temp, tp_temp, dr_temp] = deal(zeros(1, stim_comb));

for i = 1:stim_comb
    sp_temp(i) = datarun.stimulus.trials(i).SPATIAL_PERIOD;
    tp_temp(i) = datarun.stimulus.trials(i).TEMPORAL_PERIOD;
    dr_temp(i) = datarun.stimulus.trials(i).DIRECTION;
end

i_stim = zeros(length(sp), length(tp), length(dr));
index = zeros(length(sp), length(tp), length(dr), repeat);

for s = 1:length(sp)
    for t = 1:length(tp)
        for d = 1:length(dr)
            [~, i_s] = find(sp_temp == sp(s));
            [~, i_t] = find(tp_temp == tp(t));
            [~, i_d] = find(dr_temp == dr(d));
            i_sp = intersect(i_s, i_t);
            i = intersect(i_sp, i_d);
            i_stim(s, t, d) = i;
            index(s, t, d, :) = index_temp(i, :);

        end
    end
end


    

