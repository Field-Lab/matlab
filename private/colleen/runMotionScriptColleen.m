%runMotionScriptMalcom()




data_set = {'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1''2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1''2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1''2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1'};

data_run = [12;12;12;12;15;15;15;15;18;18;18;18;16;16;16;16;13;13;13;13;19;19;19;19]
config_num = [14;14;14;14;3;3;3;3;1;1;1;1;4;4;4;4;1;1;1;1;1;1;1;1];
cell_type = {'Off parasol'; 'On midget'; 'Off midget'; 'On parasol';'Off parasol'; 'On midget'; 'Off midget'; 'On parasol';'Off parasol'; 'On midget'; 'Off midget'; 'On parasol';'Off parasol'; 'On midget'; 'Off midget'; 'On parasol';'Off parasol'; 'On midget'; 'Off midget'; 'On parasol';'Off parasol'; 'On midget'; 'Off midget'; 'On parasol'};
trial_start_estimate = [118;118;118;118;118;118;118;118;118;118;118;118;220;220;220;220;220;220;220;220;220;220;220;220];


% data_set = {'2007-03-27-1';'2007-03-27-1';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4'};
% data_run = [17;17;2;2;2;3;3;3;8;8;8];
% config_num = [4;4;9;9;9;23;23;23;3;3;3];
% cell_type = {'On midget'; 'On parasol'; 'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol'};
% trial_start_estimate = [12;12;96;96;96;200;200;200;12;12;12];

for i = 24:length(config_num)
motion_script_colleen_asFunction(data_set{i}, data_run(i), config_num(i), cell_type{i}, trial_start_estimate(i))
fprintf('Trial Number %d done\n',i )
end

gmail('crhoades227@gmail.com', 'Bertha done')
