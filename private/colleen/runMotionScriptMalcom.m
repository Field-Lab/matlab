%runMotionScriptMalcom()


data_set = {'2007-03-27-1';'2007-03-27-1';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4'};
data_run = [17;17;2;2;2;3;3;3;5;8;8;8];
config_num = [4;4;9;9;9;23;23;23;3;3;3];
cell_type = {'On midget'; 'On parasol'; 'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol'};
trial_start_estimate = [12;12;96;96;96;200;200;200;12;12;12];
for i = 1:length(config_num)
motion_script_malcolm_asFunction(data_set{i}, data_run(i), config_num(i), cell_type{i}, trial_start_estimate(i))
fprintf('Trial Number %d done\n',i )
end