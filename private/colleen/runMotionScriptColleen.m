%runMotionScriptMalcom()



data_set = {'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';;'2007-08-24-4';'2007-08-24-4';'2007-08-24-4'};
data_run = [9;9;9;10;10;10;11;11;11;6;6;6;7;7;7;4;4;4;5;5;5];
config_num = [1;1;1;2;2;2;1;1;1;3;3;3;3;3;3;1;1;1;4;4;4];
cell_type = {'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol'; 'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol'};
trial_start_estimate = [12;12;12;50;50;50;50;50;50;96;96;96;96;96;96;195;195;195;195;195;195];

% data_set = {'2007-03-27-1';'2007-03-27-1';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4'};
% data_run = [17;17;2;2;2;3;3;3;8;8;8];
% config_num = [4;4;9;9;9;23;23;23;3;3;3];
% cell_type = {'On midget'; 'On parasol'; 'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol'};
% trial_start_estimate = [12;12;96;96;96;200;200;200;12;12;12];

for i = 1:length(config_num)
motion_script_colleen_asFunction(data_set{i}, data_run(i), config_num(i), cell_type{i}, trial_start_estimate(i))
fprintf('Trial Number %d done\n',i )
end

gmail('crhoades227@gmail.com', 'Bertha done')
