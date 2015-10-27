function index = grp_stim_mb(datarun)

% index = grp_stim_mb(datarun)
% group stimulus trials by stimulus combination, return the index
% index: M bar width x N speed x R directions x S repeats
%
% xyao
% 2014-10-24

list = datarun.stimulus.trial_list;
repeat = datarun.stimulus.repetitions;
stim_comb = length(datarun.stimulus.combinations);

bw = datarun.stimulus.params.BAR_WIDTH;
dt = datarun.stimulus.params.DELTA;
dr = datarun.stimulus.params.DIRECTION;



index_temp = zeros(stim_comb, repeat);

for i = 1:stim_comb
    index_temp(i, :) = find(list == i);
end

[sp_temp, tp_temp, dr_temp] = deal(zeros(1, stim_comb));

for i = 1:stim_comb
    sp_temp(i) = datarun.stimulus.trials(i).BAR_WIDTH;
    tp_temp(i) = datarun.stimulus.trials(i).DELTA;
    dr_temp(i) = datarun.stimulus.trials(i).DIRECTION;
end

i_stim = zeros(length(bw), length(dt), length(dr));
index = zeros(length(bw), length(dt), length(dr), repeat);

for s = 1:length(bw)
    for t = 1:length(dt)
        for d = 1:length(dr)
            [~, i_s] = find(sp_temp == bw(s));
            [~, i_t] = find(tp_temp == dt(t));
            [~, i_d] = find(dr_temp == dr(d));
            i_sp = intersect(i_s, i_t);
            i = intersect(i_sp, i_d);
            i_stim(s, t, d) = i;
            index(s, t, d, :) = index_temp(i, :);

        end
    end
end


    

