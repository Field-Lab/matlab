%runMotionScriptMalcom()




% data_set = {'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4'};
% data_set = {'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4'};
data_set = {'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1'; '2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4'};
% data_set = {'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4'};
data_run = [12;12;12;12;13;13;13;13;15;15;15;15;16;16;16;16;18;18;18;18;19;19;19;19;2;2;2;3;3;3;4;4;4;5;5;5;6;6;6;7;7;7;10;10;10;11;11;11];
% data_run = [16;16;16;16];
% config_num = [14;14;14;14;1;1;1;1;3;3;3;3;3;3;3;3;4;4;4;4;4;4;4;4;1;1;1;1];
config_num = [16;16;16;16;10;10;10;10;2;2;2;2;3;3;3;3;3;3;3;3;2;2;2;2;19;19;19;17;17;17;2;2;2;3;3;3;1;1;1;4;4;4;3;3;3;2;2;2];
% cell_type = {'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';};
% cell_type = {'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol'};

cell_type = {'Off parasol'; 'On midget'; 'On parasol';'Off midget'; 'Off parasol'; 'On midget'; 'On parasol';'Off midget';'Off parasol'; 'On midget'; 'On parasol';'Off midget';  'Off parasol'; 'On midget'; 'On parasol';'Off midget'; 'Off parasol'; 'On midget'; 'On parasol';'Off midget'; 'Off parasol'; 'On midget'; 'On parasol';'Off midget';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol'};
% trial_start_estimate = [230;230;230;230];

% data_set = {'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4'};
% data_run = [17;17;2;2;2;3;3;3;8;8;8];
% config_num = [4;4;9;9;9;23;23;23;3;3;3];
% cell_type = {'On midget'; 'On parasol'; 'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol'};
% trial_start_estimate = [12;12;96;96;96;200;200;200;12;12;12];

velocity = 12*[8;8;8;8;16;16;16;16;8;8;8;8;16;16;16;16;8;8;8;8;16;16;16;16;8;8;8;16;16;16;16;16;16;16;16;16;8;8;8;8;8;8;4;4;4;4;4;4];

for i =2:length(config_num)
motion_script_colleen_asFunction(data_set{i}, data_run(i), config_num(i), cell_type{i}, velocity(i))
fprintf('Trial Number %d done\n',i )
end

gmail('crhoades227@gmail.com', 'Bertha done')