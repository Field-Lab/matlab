%% Formatting units for easy processing (description look in 'S:\data\alexandra\for_Thomas\for_Al\spike_data_Al\ReadMe.txt')

cd('/mnt/muench_data/user/alexandra/scripts')
% check for the time unit (in my settings are seconds.)
date='20120405_1';
a=dir(['/mnt/muench_data/data/alexandra/MEA_data/',date,'/units/*.mat']);
for i=1:length(a)
    load(['/mnt/muench_data/data/alexandra/MEA_data/',date,'/units/',a(i).name])
    if strcmp(unit{1}(5,2),'1000')
        display(i)
        unit{1}(5,2)={'1'};
        save(['/mnt/muench_data/data/alexandra/MEA_data/',date,'/units/',a(i).name])
    end
end


%% Checker Board Flicker spike-info (not shuffling)
cd('/mnt/muench_data/user/alexandra/scripts')
date='20120902_1';

mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date];
hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

if ~exist([mainpath '/easy_formatted_units/'],'dir')    
    mkdir([mainpath, '/easy_formatted_units']);
end
path2save=[mainpath '/easy_formatted_units/'];

% make list of files to put into easy formatted unit (choosing by the name of the stimulus in heka files)
file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,'check_seq1%bg_1%_', 'once'))
        file_list=[file_list i];
    end
end

timings=dir([mainpath,'/timing/*mat']);

%get list of units to be reformatted
basic_format_units_list=dir([mainpath, '/units/*.mat']);

clear spike_info
for cnt=1:length(basic_format_units_list)
    
    basic_format_units_list(cnt).name
    
    load([mainpath,'/units/',basic_format_units_list(cnt).name]);
    for i=1:length(file_list)
        protocol=read_header_field_heka(hekapath, heka(file_list(i)).name, 'Stimulus Protocol');
        load([mainpath,'/timing/',timings(i).name])
        FlipTimeStamps=[(FlipTimeStamps-FlipTimeStamps(1))'*1000+protocol(find((protocol(:,3)==9))-1,1); protocol(find((protocol(:,3)==9)):end,1)]; % in ms; time of the stimulation begin added, time of the end of stimulation and background presentation appended
        
        spike_info.spike_times{i}=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        spike_info.flip_times{i}=[round(FlipTimeStamps) [-1; (1:length(FlipTimeStamps)-3)'; -1; -1]];
        spike_info.name_info{i}=[basic_format_units_list(cnt).name(1:end-4),'___',heka(file_list(i)).name(13:end-5)];
    end
    if size(spike_info.spike_times,1)<size(spike_info.spike_times,2)
        spike_info.spike_times=spike_info.spike_times';
    end
    if size(spike_info.flip_times,1)<size(spike_info.flip_times,2)
        spike_info.flip_times=spike_info.flip_times';
    end
    if size(spike_info.name_info,1)<size(spike_info.name_info,2)
        spike_info.name_info=spike_info.name_info';
    end
    
    save([path2save,basic_format_units_list(cnt).name(1:end-4),'___spike_info_CBFlicker'],'spike_info')
end



%% Checker Board Flicker spike-info (single cone)
cd('/mnt/muench_data/user/alexandra/scripts')
date='20121024';

mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date];
hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

if ~exist([mainpath '/easy_formatted_units/'],'dir')    
    mkdir([mainpath, '/easy_formatted_units']);
end
path2save=[mainpath '/easy_formatted_units/'];

% make list of files to put into easy formatted unit (choosing by the name of the stimulus in heka files)
file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,'rng', 'once'))
        file_list=[file_list i];
    end
end

timings=dir([mainpath,'/timing_sg/*mat']);

%get list of units to be reformatted
basic_format_units_list=dir([mainpath, '/units/*.mat']);

clear spike_info
for cnt=1:length(basic_format_units_list)
    
    basic_format_units_list(cnt).name
    
    load([mainpath,'/units/',basic_format_units_list(cnt).name]);
    for i=1:length(file_list)
        protocol=read_header_field_heka(hekapath, heka(file_list(i)).name, 'Stimulus Protocol');
        load([mainpath,'/timing_sg/',timings(i).name])
        FlipTimeStamps=[(FlipTimeStamps-FlipTimeStamps(1))'*1000+protocol(find((protocol(:,3)==6))-1,1); protocol(find((protocol(:,3)==9)):end,1)]; % in ms; time of the stimulation begin added, time of the end of stimulation and background presentation appended
        
        spike_info.spike_times{i}=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        spike_info.flip_times{i}=[round(FlipTimeStamps) [-1; (1:length(FlipTimeStamps)-3)'; -1; -1]];
        spike_info.name_info{i}=[basic_format_units_list(cnt).name(1:end-4),'___',heka(file_list(i)).name(13:end-5)];
    end
    if size(spike_info.spike_times,1)<size(spike_info.spike_times,2)
        spike_info.spike_times=spike_info.spike_times';
    end
    if size(spike_info.flip_times,1)<size(spike_info.flip_times,2)
        spike_info.flip_times=spike_info.flip_times';
    end
    if size(spike_info.name_info,1)<size(spike_info.name_info,2)
        spike_info.name_info=spike_info.name_info';
    end
    
    save([path2save,basic_format_units_list(cnt).name(1:end-4),'___spike_info_SingleCone'],'spike_info')
end



%% Checker Board Flicker spike-info (shuffling)
date='20120927';

mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date];
hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

if ~exist([mainpath 'easy_formatted_units/'],'dir')    
    mkdir([mainpath, 'easy_formatted_units']);
end
path2save=[mainpath '/easy_formatted_units/'];

% make list of files to put into easy formatted unit (choosing by the name of the stimulus in heka files)
file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,'check_seq1%bg_1_', 'once'))
        file_list=[file_list i];
    end
end

timings=dir([mainpath,'/timing/*mat']);

%get list of units to be reformatted
basic_format_units_list=dir([mainpath, '/units/*.mat']);

clear spike_info
for cnt=1:length(basic_format_units_list)
    
    basic_format_units_list(cnt).name
    
    load([mainpath,'/units/',basic_format_units_list(cnt).name]);
    for i=1:length(file_list)
        protocol=read_header_field_heka(hekapath, heka(file_list(i)).name, 'Stimulus Protocol');
        load([mainpath,'/timing/',timings(i).name])
        FlipTimeStamps=[(FlipTimeStamps-FlipTimeStamps(1))'*1000+protocol(find((protocol(:,3)==11))-1,1); protocol(find((protocol(:,3)==11)):end,1)]; % in ms; time of the stimulation begin added, time of the end of stimulation and background presentation appended
        
        spike_info.spike_times{i}=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        spike_info.flip_times{i}=[round(FlipTimeStamps) [-1; (1:length(FlipTimeStamps)-3)'; -1; -1]];
        spike_info.name_info{i}=[basic_format_units_list(cnt).name(1:end-4),'___',heka(file_list(i)).name(13:end-5)];
    end
    if size(spike_info.spike_times,1)<size(spike_info.spike_times,2)
        spike_info.spike_times=spike_info.spike_times';
    end
    if size(spike_info.flip_times,1)<size(spike_info.flip_times,2)
        spike_info.flip_times=spike_info.flip_times';
    end
    if size(spike_info.name_info,1)<size(spike_info.name_info,2)
        spike_info.name_info=spike_info.name_info';
    end
    
    save([mainpath,'easy_formatted_units/',basic_format_units_list(cnt).name(1:end-4),'___spike_info_CBFlicker_Shuffle'],'spike_info')
end


%% Full Field Flicker spike-info
cd('/mnt/muench_data/user/alexandra/scripts')
date='20130226';

mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

if ~exist([mainpath 'easy_formatted_units/'],'dir')    
    mkdir([mainpath, 'easy_formatted_units']);
end
path2save=[mainpath 'easy_formatted_units/'];

% make list of files to put into easy formatted unit (choosing by the name of the stimulus in heka files)
file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,'Cseq', 'once'))
        file_list=[file_list i];
    end
end

%get list of units to be reformatted
basic_format_units_list=dir([mainpath, 'units/*.mat']);

% get protocol
protocol=read_header_field_heka(hekapath, heka(file_list(1)).name, 'Stimulus Protocol');
protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
protocol=round(protocol(2:end,[1,4]));
protocol_accum=zeros(size(protocol,1),2,length(file_list));
for i=1:length(file_list)
        protocol=read_header_field_heka(hekapath, heka(file_list(i)).name, 'Stimulus Protocol');
        protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
        protocol=round(protocol(2:end,[1,4]));
        protocol_accum(:,:,i)=protocol;
end

for cnt=1:length(basic_format_units_list)
    basic_format_units_list(cnt).name
    load([mainpath,'units/',basic_format_units_list(cnt).name]);
    for i=1:length(file_list)
        spike_info.spike_times{i}=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms                
        spike_info.flip_times{i}=protocol_accum(:,:,i);
        spike_info.name_info{i}=[basic_format_units_list(cnt).name(1:end-4),'___',heka(file_list(i)).name(13:end-5)];
    end
    if size(spike_info.spike_times,1)<size(spike_info.spike_times,2)
        spike_info.spike_times=spike_info.spike_times';
    end
    if size(spike_info.flip_times,1)<size(spike_info.flip_times,2)
        spike_info.flip_times=spike_info.flip_times';
    end
    if size(spike_info.name_info,1)<size(spike_info.name_info,2)
        spike_info.name_info=spike_info.name_info';
    end
    save([mainpath,'easy_formatted_units/',basic_format_units_list(cnt).name(1:end-4),'___spike_info_FFFlicker'],'spike_info')
end


%% Full Field Flash spike-info
clear
cd('/mnt/muench_data/user/alexandra/scripts')
date='20130228_1';

mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

if ~exist([mainpath 'easy_formatted_units/'],'dir')    
    mkdir([mainpath, 'easy_formatted_units']);
end
path2save=[mainpath 'easy_formatted_units/'];

file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,'quick', 'once')) 
        file_list=[file_list i];
    end
end

%get list of units to be reformatted
basic_format_units_list=dir([mainpath, 'units/*.mat']);

for cnt=1:length(basic_format_units_list)
    basic_format_units_list(cnt).name
    load([mainpath,'units/',basic_format_units_list(cnt).name]);
    for i=1:length(file_list)
        protocol=read_header_field_heka(hekapath, heka(file_list(i)).name, 'Stimulus Protocol');
        spike_info.spike_times{i}=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        protocol(protocol(:,2)==1,2)=50;
        protocol(protocol(:,2)==0,2)=10;
        spike_info.flip_times{i}=round(protocol(2:end,[1,2]));
        spike_info.name_info{i}=[basic_format_units_list(cnt).name(1:end-4),'___',heka(file_list(i)).name(13:end-5)];
    end
    if size(spike_info.spike_times,1)<size(spike_info.spike_times,2)
        spike_info.spike_times=spike_info.spike_times';
    end
    if size(spike_info.flip_times,1)<size(spike_info.flip_times,2)
        spike_info.flip_times=spike_info.flip_times';
    end
    if size(spike_info.name_info,1)<size(spike_info.name_info,2)
        spike_info.name_info=spike_info.name_info';
    end

    save([mainpath,'easy_formatted_units/',basic_format_units_list(cnt).name(1:end-4),'___spike_info_FFFlash'],'spike_info')
end

%% Moment Flash spike-info
clear
cd('/mnt/muench_data/user/alexandra/scripts')
date='20130226';

mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

if ~exist([mainpath 'easy_formatted_units/'],'dir')    
    mkdir([mainpath, 'easy_formatted_units']);
end
path2save=[mainpath 'easy_formatted_units/'];

file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,'moment', 'once')) 
        file_list=[file_list i];
    end
end

%get list of units to be reformatted
basic_format_units_list=dir([mainpath, 'units/*.mat']);

for cnt=1:length(basic_format_units_list)
    basic_format_units_list(cnt).name
    load([mainpath,'units/',basic_format_units_list(cnt).name]);
    for i=1:length(file_list)
        protocol=read_header_field_heka(hekapath, heka(file_list(i)).name, 'Stimulus Protocol');
        spike_info.spike_times{i}=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        protocol(protocol(:,2)==1,2)=50;
        protocol(protocol(:,2)==0,2)=10;
        spike_info.flip_times{i}=round(protocol(2:end,[1,2]));
        spike_info.name_info{i}=[basic_format_units_list(cnt).name(1:end-4),'___',heka(file_list(i)).name(13:end-5)];
    end
    if size(spike_info.spike_times,1)<size(spike_info.spike_times,2)
        spike_info.spike_times=spike_info.spike_times';
    end
    if size(spike_info.flip_times,1)<size(spike_info.flip_times,2)
        spike_info.flip_times=spike_info.flip_times';
    end
    if size(spike_info.name_info,1)<size(spike_info.name_info,2)
        spike_info.name_info=spike_info.name_info';
    end

    save([mainpath,'easy_formatted_units/',basic_format_units_list(cnt).name(1:end-4),'___spike_info_Moment'],'spike_info')
end


%% ND spike-info
clear
cd('/mnt/muench_data/user/alexandra/scripts')
date='20130228_1';

mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

if ~exist([mainpath 'easy_formatted_units/'],'dir')    
    mkdir([mainpath, 'easy_formatted_units']);
end
path2save=[mainpath 'easy_formatted_units/'];

file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,'nd', 'once')) 
        file_list=[file_list i];
    end
end

%get list of units to be reformatted
basic_format_units_list=dir([mainpath, 'units/*.mat']);

for cnt=1:length(basic_format_units_list)
    basic_format_units_list(cnt).name
    load([mainpath,'units/',basic_format_units_list(cnt).name]);
    for i=1:length(file_list)
        protocol=read_header_field_heka(hekapath, heka(file_list(i)).name, 'Stimulus Protocol');
        spike_info.spike_times{i}=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        spike_info.flip_times{i}=round(protocol(2:end,[1,4]));
        spike_info.name_info{i}=[basic_format_units_list(cnt).name(1:end-4),'___',heka(file_list(i)).name(13:end-5)];
    end
    if size(spike_info.spike_times,1)<size(spike_info.spike_times,2)
        spike_info.spike_times=spike_info.spike_times';
    end
    if size(spike_info.flip_times,1)<size(spike_info.flip_times,2)
        spike_info.flip_times=spike_info.flip_times';
    end
    if size(spike_info.name_info,1)<size(spike_info.name_info,2)
        spike_info.name_info=spike_info.name_info';
    end

    save([mainpath,'easy_formatted_units/',basic_format_units_list(cnt).name(1:end-4),'___spike_info_ND'],'spike_info')
end




%% Chirp spike-info
clear
cd('/mnt/muench_data/user/alexandra/scripts')
date='20130228_1';

mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

if ~exist([mainpath 'easy_formatted_units/'],'dir')    
    mkdir([mainpath, 'easy_formatted_units']);
end
path2save=[mainpath 'easy_formatted_units/'];

file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,'chirp', 'once')) 
        file_list=[file_list i];
    end
end

%get list of units to be reformatted
basic_format_units_list=dir([mainpath, 'units/*.mat']);

protocol=read_header_field_heka(hekapath, heka(file_list(1)).name, 'Stimulus Protocol');
protocol(protocol(:,2)>0,2)=protocol(protocol(:,2)>0,2)+20;
protocol=round([protocol(2:end,1) [-1; protocol(3:end-1,4); -1]]);
protocol_accum=zeros(size(protocol,1),2,length(file_list));
for i=1:length(file_list)
    
    protocol=read_header_field_heka(hekapath, heka(file_list(i)).name, 'Stimulus Protocol');
    protocol(protocol(:,2)>0,2)=protocol(protocol(:,2)>0,2)+20;
    protocol=round([protocol(2:end,1) [-1; protocol(3:end-1,4); -1]]);
    protocol_accum(:,:,i)=protocol;
end
        

for cnt=1:length(basic_format_units_list)
    
    basic_format_units_list(cnt).name
    load([mainpath,'units/',basic_format_units_list(cnt).name]);
    for i=1:length(file_list)
        spike_info.spike_times{i}=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        spike_info.flip_times{i}=protocol_accum(:,:,i);
        spike_info.name_info{i}=[basic_format_units_list(cnt).name(1:end-4),'___',heka(file_list(i)).name(13:end-5)];
    end
    if size(spike_info.spike_times,1)<size(spike_info.spike_times,2)
        spike_info.spike_times=spike_info.spike_times';
    end
    if size(spike_info.flip_times,1)<size(spike_info.flip_times,2)
        spike_info.flip_times=spike_info.flip_times';
    end
    if size(spike_info.name_info,1)<size(spike_info.name_info,2)
        spike_info.name_info=spike_info.name_info';
    end
    
    save([mainpath,'easy_formatted_units/',basic_format_units_list(cnt).name(1:end-4),'___spike_info_chirp'],'spike_info')
end



%% Full Field Flicker HL5s spike-info
clear
cd('/mnt/muench_data/user/alexandra/scripts')
date='20130228_1';
codeWord='sine';

mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

if ~exist([mainpath 'easy_formatted_units/'],'dir')    
    mkdir([mainpath, 'easy_formatted_units']);
end
path2save=[mainpath 'easy_formatted_units/'];

% make list of files to put into easy formatted unit (choosing by the name of the stimulus in heka files)
file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,codeWord, 'once'))
        file_list=[file_list i];
    end
end

timingPath=[mainpath 'timing/'];
timings=dir([timingPath,'*.mat']);

%get list of units to be reformatted
basic_format_units_list=dir([mainpath, 'units/*.mat']);

% get protocol
protocol=read_header_field_heka(hekapath, heka(file_list(1)).name, 'Stimulus Protocol');
protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
protocol=round(protocol(2:end,[1,4]));
protocol_accum=zeros(size(protocol,1),2,length(file_list));
abused_accum=zeros(size(protocol,1),length(file_list));
for i=1:length(file_list)
        protocol=read_header_field_heka(hekapath, heka(file_list(i)).name, 'Stimulus Protocol');
        protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
        protocol=round(protocol(2:end,[1,4]));
        protocol_accum(:,:,i)=protocol;
        load([timingPath, timings(file_list(i)).name]);
        abused=abused(2:end,5)-abused(2,5);
        abused_accum(:,i)=abused(3:end)*1000;
end


for cnt=1:length(basic_format_units_list)
    basic_format_units_list(cnt).name
    load([mainpath,'units/',basic_format_units_list(cnt).name]);
    for i=1:length(file_list)
        spike_info.spike_times{i}=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms                
        spike_info.flip_times{i}=[protocol_accum(:,:,i) abused_accum(:,i)];
        spike_info.name_info{i}=[basic_format_units_list(cnt).name(1:end-4),'___',heka(file_list(i)).name(13:end-5)];
    end
    if size(spike_info.spike_times,1)<size(spike_info.spike_times,2)
        spike_info.spike_times=spike_info.spike_times';
    end
    if size(spike_info.flip_times,1)<size(spike_info.flip_times,2)
        spike_info.flip_times=spike_info.flip_times';
    end
    if size(spike_info.name_info,1)<size(spike_info.name_info,2)
        spike_info.name_info=spike_info.name_info';
    end
    save([mainpath,'easy_formatted_units/',basic_format_units_list(cnt).name(1:end-4),'___spike_info_',codeWord],'spike_info')
end
