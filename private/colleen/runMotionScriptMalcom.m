%runMotionScriptMalcom()


data_set = {'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1''2007-03-27-1''2007-03-27-1';'2007-03-27-1';'2007-08-24-4';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1''2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1'};
data_run = [12;12;12;12;16;16;16;16;17;17;17;17;2;2;2;3;3;3;5;5;5;7;7;7;8;8;8];
config_num = [14;14;14;14;4;4;4;4;4;4;4;4;9;9;9;23;23;23;4;4;4;3;3;3;3;3;3];
cell_type = {'Off midget'; 'Off parasol'; 'On midget'; 'On parasol'; 'Off midget'; 'Off parasol'; 'On midget'; 'On parasol'; 'Off midget'; 'Off parasol'; 'On midget'; 'On parasol'; 'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol'};

for i = 1:length(config_num)
motion_script_malcolm_asFunction(data_set{i}, data_run(i), config_num(i), cell_type{i})
end