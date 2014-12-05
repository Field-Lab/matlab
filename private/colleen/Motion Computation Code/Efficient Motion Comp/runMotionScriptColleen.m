% runMotionScriptColleen
% Set parameters for motion_script_colleen_asFunction and be able to
% compute estimates for a lot of runs in a single program


data_set = {'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1'; '2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4'};
data_run = [12;12;12;12;13;13;13;13;15;15;15;15;16;16;16;16;18;18;18;18;19;19;19;19;2;2;2;3;3;3;4;4;4;5;5;5;6;6;6;7;7;7;10;10;10;11;11;11];

% for which stimulus to run
config_num = [14;14;14;14;1;1;1;1;3;3;3;3;4;4;4;4;1;1;1;1;1;1;1;1;9;9;9;23;23;23;1;1;1;4;4;4;3;3;3;3;3;3;2;2;2;1;1;1];

cell_type = {'Off parasol'; 'On midget'; 'On parasol';'Off midget'; 'Off parasol'; 'On midget'; 'On parasol';'Off midget';'Off parasol'; 'On midget'; 'On parasol';'Off midget';  'Off parasol'; 'On midget'; 'On parasol';'Off midget'; 'Off parasol'; 'On midget'; 'On parasol';'Off midget'; 'Off parasol'; 'On midget'; 'On parasol';'Off midget';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol'};

for i =1:length(config_num)
    motion_script_colleen_asFunction(data_set{i}, data_run(i), config_num(i), cell_type{i}, velocity(i))
    fprintf('Trial Number %d done\n',i )
end

% Send email when complete
gmail('crhoades227@gmail.com', 'my computer done')
