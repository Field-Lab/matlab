function get_LFP_linux_all(ch,date)

pathway=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/MEA_bin/channels/channel',int2str(ch),'/'];
file_list=dir([pathway,'*.mat']);

path2save=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LFP_all/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

if ~exist([path2save,'channel',int2str(ch)],'file')
        mkdir([path2save,'channel',int2str(ch)]);
end
path2save=[path2save,'channel',int2str(ch),'/'];

already=dir([path2save,'*.mat']);

% LOW PASS - for LFP
sf=25000;%sampling frequency
passL=100;orderL=2;
fNormL = passL/(sf/2);
[aL,bL] = butter(orderL, fNormL, 'low');


for fileID=1:length(file_list)
    if ~strcmp(already(fileID).name,file_list(fileID).name)
        load([pathway, file_list(fileID).name]);        
        lfp = filtfilt(aL, bL, data');
        lfp=lfp(1:25:end);
        save([path2save,file_list(fileID).name],'lfp');
    end
end




