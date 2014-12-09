% runMotionScriptColleen
% Set parameters for motion_script_colleen_asFunction and be able to
% compute estimates for a lot of runs in a single program


data_set = {'2007-03-27-1';'2007-03-27-1';'2007-03-27-1'; '2007-03-27-1';'2007-03-27-1';'2007-03-27-1';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4';'2007-08-24-4'};
data_run = [12;13;15;16;18;19;2;3;4;5;6;7;10;11];

% for which stimulus to run
config_num = [14;1;3;4;1;1;9;23;1;4;3;3;2;1];

cell_type = {'Off parasol'; 'On midget'; 'On parasol';'Off midget'; 'Off parasol'; 'On midget'; 'On parasol';'Off midget';'Off parasol'; 'On midget'; 'On parasol';'Off midget';  'Off parasol'; 'On midget'; 'On parasol';'Off midget'; 'Off parasol'; 'On midget'; 'On parasol';'Off midget'; 'Off parasol'; 'On midget'; 'On parasol';'Off midget';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol';'Off parasol'; 'On midget'; 'On parasol'};
velocity = [96;192;96;192;96;192;96;192;192;192;96;96;48;48];
for i =2:length(config_num)
    motion_script_pooled_asFunction(data_set{i}, data_run(i), config_num(i), velocity(i))
    fprintf('Trial Number %d done\n',i )
end

% Send email when complete
gmail('crhoades227@gmail.com', 'bertha done')
