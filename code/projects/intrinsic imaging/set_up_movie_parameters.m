% SET UP PARAMETERS


switch 8
    case 1 % 2010-02-26-0/data000

        cell_id = 6811;

        datarun = load_data('2010-02-26-0/data000/data000');


        % load triggers and spikes
        datarun = load_neurons(datarun);
        spikes = datarun.spikes{get_cell_indices(datarun,cell_id)};


        % image information
        image_spec = struct;
%        image_spec.path = '/snle/data/2010-02-26-0/spot/sequences/';
        image_spec.path = '/snle/home/peterli/Desktop/2010-02-26-0/';
        image_spec.data = struct;


        % trigger times

        seq{1}.triggers = datarun.triggers(5:4004);
        seq{1}.name = 's1.tif';

        seq{2}.triggers = datarun.triggers(4006:8005);
        seq{2}.name = 's2.tif';

        seq{3}.triggers = datarun.triggers(8007:12006);
        seq{3}.name = 's3.tif';

        seq{4}.triggers = datarun.triggers(12008:16007);
        seq{4}.name = 's4.tif';


        seq_data = cell(length(seq));
        for ss=1:length(seq)
            seq_data{ss}=struct;
            for ii=1:length(seq{ss}.triggers)
                seq_data{ss}(ii).name = seq{ss}.name;
                seq_data{ss}(ii).index = ii;
                seq_data{ss}(ii).time = seq{ss}.triggers(ii);
            end
        end


        % choose sequence
        switch 1
            case 0 % everything
                image_spec.data = [seq_data{:}];
            case 1  % s1
                image_spec.data = seq_data{1};
            case 2  % s2
                image_spec.data = seq_data{2};
            case 3  % s1 subset
                image_spec.data = seq_data{1}(1:200);
        end




    case 2 % 2010-03-01-0/data002

        %cell_id = 7235;

        datarun = load_data('2010-03-01-0/data002/data002');




        % load triggers
        datarun = load_neurons(datarun);
        triggers = datarun.triggers([1:4000 4002:8001 8003:12002 12004:16003]);
        %triggers = triggers(randperm(length(triggers)));


        % get spikes
        %spikes = datarun.spikes{get_cell_indices(datarun,cell_id)};
        spikes = -0.32:8:datarun.duration;

        % image information
        image_spec = struct;
        image_spec.path = '/snle/data/2010-03-01-0/spot/data002/';
        image_spec.data = struct;

        [image_spec.data(1:4000).name] = deal('s1.tif');
        for ii=1:4000;
            image_spec.data(ii).index = ii;
            image_spec.data(ii).time = triggers(ii);
        end

        [image_spec.data(4001:8000).name] = deal('s2.tif');
        for ii=4001:8000;
            image_spec.data(ii).index = ii-4000;
            image_spec.data(ii).time = triggers(ii);
        end

        [image_spec.data(8001:12000).name] = deal('s3.tif');
        for ii=8001:12000;
            image_spec.data(ii).index = ii-8000;
            image_spec.data(ii).time = triggers(ii);
        end

        [image_spec.data(12001:16000).name] = deal('s4.tif');
        for ii=12001:16000;
            image_spec.data(ii).index = ii-12000;
            image_spec.data(ii).time = triggers(ii);
        end





        % choose sequence
        switch 2
            case 0 % everything
            case 1  % s1
                image_spec.data = image_spec.data(1:4000);
            case 2  % s2
                image_spec.data = image_spec.data(4001:8000);
            case 3  % s3
                image_spec.data = image_spec.data(8001:12000);
            case 4  % s4
                image_spec.data = image_spec.data(12001:16000);

            case 5  % s1 subset
                image_spec.data = image_spec.data(1:100);
        end




    case 3 % 2010-03-01-0/data004

        %cell_id = 6485;
        cell_id = 6496;

        datarun = load_data('2010-03-01-0/data004/data004');


        % load triggers
        datarun = load_neurons(datarun);
        triggers = datarun.triggers([   2:4001   4005:8004   8005:12004   12008:16007   16009:20008   20009:24008   ]);

        % get spikes
        spikes = datarun.spikes{get_cell_indices(datarun,cell_id)};

        % image information
        image_spec = struct;
        image_spec.path = '/snle/data/2010-03-01-0/spot/data004/';
        image_spec.data = struct;

        [image_spec.data(1:4000).name] = deal('s1.tif');
        for ii=1:4000;
            image_spec.data(ii).index = ii;
            image_spec.data(ii).time = triggers(ii);
        end

        [image_spec.data(4001:8000).name] = deal('s2.tif');
        for ii=4001:8000;
            image_spec.data(ii).index = ii-4000;
            image_spec.data(ii).time = triggers(ii);
        end

        [image_spec.data(8001:12000).name] = deal('s3.tif');
        for ii=8001:12000;
            image_spec.data(ii).index = ii-8000;
            image_spec.data(ii).time = triggers(ii);
        end

        [image_spec.data(12001:16000).name] = deal('s4.tif');
        for ii=12001:16000;
            image_spec.data(ii).index = ii-12000;
            image_spec.data(ii).time = triggers(ii);
        end

        [image_spec.data(16001:20000).name] = deal('s5.tif');
        for ii=16001:20000;
            image_spec.data(ii).index = ii-16000;
            image_spec.data(ii).time = triggers(ii);
        end

        [image_spec.data(20001:24000).name] = deal('s6.tif');
        for ii=20001:24000;
            image_spec.data(ii).index = ii-20000;
            image_spec.data(ii).time = triggers(ii);
        end





        % choose sequence
        switch 2
            case 0 % everything
            case 1  % s1-s3 (focused on electrodes)
                image_spec.data = image_spec.data(1:12000);
            case 2  % s4-s6 (focused 5µm above electrodes)
                image_spec.data = image_spec.data(12001:24000);
            case 3  % subset of s1
                image_spec.data = image_spec.data(1:2000);
            case 4  % subset of s1 with apparent effect for cell 6485
                image_spec.data = image_spec.data(1:400);
        end




    case 4 % 2010-03-01-0/data005

        %cell_id = 5133; % big on 344
        cell_id = 5331; % big on 33
        %cell_id = 603; % big on 41
        %cell_id = 604; % big on 41
        %cell_id = 606; % medium on 41

        datarun = load_data('2010-03-01-0/data005/data005');


        % load triggers, spikes
        datarun = load_neurons(datarun);

        % get spikes
        spikes = datarun.spikes{get_cell_indices(datarun,cell_id)};

        % image information
        image_spec = struct;
        image_spec.path = '/snle/data/2010-03-01-0/spot/data005/';
        image_spec.data = struct;


        seq{1}.triggers = datarun.triggers(1:4000);
        seq{1}.name = 's1.tif';

        seq{2}.triggers = datarun.triggers(4001:8000);
        seq{2}.name = 's2.tif';

        seq{3}.triggers = datarun.triggers(8001:12000);
        seq{3}.name = 's3.tif';

        seq{4}.triggers = datarun.triggers(12005:16004);
        seq{4}.name = 's4.tif';

        seq{5}.triggers = datarun.triggers(16009:20008);
        seq{5}.name = 's5.tif';

        seq{6}.triggers = datarun.triggers(20017:22919);
        seq{6}.name = 's6.tif';

        seq{7}.triggers = datarun.triggers(22920:25822);
        seq{7}.name = 's7.tif';


        seq_data = cell(length(seq));
        for ss=1:length(seq)
            seq_data{ss}=struct;
            for ii=1:length(seq{ss}.triggers)
                seq_data{ss}(ii).name = seq{ss}.name;
                seq_data{ss}(ii).index = ii;
                seq_data{ss}(ii).time = seq{ss}.triggers(ii);
            end
        end


        % choose sequence
        switch 3
            case 0 % everything
                image_spec.data = [seq_data{:}];
            case 1  % s6
                image_spec.data = seq_data{6};
            case 2  % s7
                image_spec.data = seq_data{7};
            case 3  % s6 subset
                image_spec.data = seq_data{6}(1:200);
        end



    case 5 % 2010-03-01-0/data006

        %cell_id = 5133; % big on 344
        cell_id = 5327; % moderate on 348
        %cell_id = 363; % big on 33
        %cell_id = 603; % big on 41
        %cell_id = 604; % big on 41

        datarun = load_data('2010-03-01-0/data006/data006');


        % load triggers, spikes
        datarun = load_neurons(datarun);
        spikes = datarun.spikes{get_cell_indices(datarun,cell_id)};

        % image information
        image_spec = struct;
        image_spec.path = '/snle/data/2010-03-01-0/spot/data006/';
        image_spec.data = struct;


        seq{1}.triggers = datarun.triggers(1:4000);
        seq{1}.name = 's1.tif';

        seq{2}.triggers = datarun.triggers(4001:8000);
        seq{2}.name = 's2.tif';

        seq{3}.triggers = datarun.triggers(8003:12002);
        seq{3}.name = 's3.tif';

        seq{4}.triggers = datarun.triggers(12003:16002);
        seq{4}.name = 's4.tif';

        seq{5}.triggers = datarun.triggers(16005:20004);
        seq{5}.name = 's5.tif';


        seq_data = cell(length(seq));
        for ss=1:length(seq)
            seq_data{ss}=struct;
            for ii=1:length(seq{ss}.triggers)
                seq_data{ss}(ii).name = seq{ss}.name;
                seq_data{ss}(ii).index = ii;
                seq_data{ss}(ii).time = seq{ss}.triggers(ii);
            end
        end


        % choose sequence
        switch 5
            case 0 % everything
                image_spec.data = [seq_data{:}];
            case 1  % s1
                image_spec.data = seq_data{1};
            case 2  % s2
                image_spec.data = seq_data{2};
            case 3  % s3
                image_spec.data = seq_data{3};
            case 4  % s4
                image_spec.data = seq_data{4};
            case 5  % s5
                image_spec.data = seq_data{5};
            case 6  % s5 subset
                image_spec.data = seq_data{5}(1:200);
        end

        
        
        

    case 6 % 2010-04-21-0/data003

        %cell_id = 7262; % near 476
        cell_id = 5506; % big on 368

        datarun = load_data('2010-04-21-0/data003/data003');


        % load triggers and spikes
        datarun = load_neurons(datarun);
        spikes = datarun.spikes{get_cell_indices(datarun,cell_id)};


        % image information
        image_spec = struct;
        image_spec.path = '/snle/data/2010-04-21-0/spot/';
        image_spec.data = struct;


        % trigger times

        % loc 1
        seq{1}.triggers = datarun.triggers(1:4000);
        seq{1}.name = 'L1/S3.tif';

        % loc 1
        seq{2}.triggers = datarun.triggers(4002:8001);
        seq{2}.name = 'L1/S4.tif';

        % loc 2
        seq{3}.triggers = datarun.triggers(8003:12002);
        seq{3}.name = 'L2/S1.tif';


        seq_data = cell(length(seq));
        for ss=1:length(seq)
            seq_data{ss}=struct;
            for ii=1:length(seq{ss}.triggers)
                seq_data{ss}(ii).name = seq{ss}.name;
                seq_data{ss}(ii).index = ii;
                seq_data{ss}(ii).time = seq{ss}.triggers(ii);
            end
        end


        % choose sequence
        switch 3
            case 0 % everything
                image_spec.data = [seq_data{:}];
            case 1  % set 1
                image_spec.data = seq_data{1};
            case 2  % set 2
                image_spec.data = seq_data{2};
            case 3  % set 2
                image_spec.data = seq_data{3};
            case 4  % subset
                image_spec.data = seq_data{3}(201:400);
        end
        
        
    case 7 % Redo case 1
        
        cell_id = 6811;
        datarun = load_data('2010-02-26-0/data000/data000');
        datarun = load_neurons(datarun);

        cell_num = get_cell_indices(datarun, cell_id);
        spikes = datarun.spikes{cell_num};
        seqs = get_camera_trigger_sequences(datarun.triggers);

        stack = struct;
        stack.basepath = '/Users/peterli/Desktop/2010-02-26-0/';
        stack.paths = repmat({'s1.tif'}, 4000, 1);
        stack.pages = num2cell((1:4000)');
        stack.triggers = seqs{1};
        
    case 8
        load intrins_2010-07-22-0_data003
        cell_id = 563;
        cell_num = get_cell_indices(datarun, cell_id);
        spikes = datarun.spikes{cell_num};
        stack = datarun.stacks{3,4};
end

clear ss ii
