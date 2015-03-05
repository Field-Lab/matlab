function get_LFP_linux(ch,date,keyWord)

pathway=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/MEA_bin/channels/channel',int2str(ch),'/'];
file_list=dir([pathway,'*',keyWord,'*.mat']);

path2save=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LFP/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end


% LOW PASS - for LFP
sf=25000;%sampling frequency
passL=100;orderL=2;
fNormL = passL/(sf/2);
[aL,bL] = butter(orderL, fNormL, 'low');

load([pathway, file_list(1).name]);
length_expected=floor(length(data)/25);

NoF=length(file_list);

lfp=zeros(length_expected,NoF);

for fileID=1:length(file_list)
    load([pathway, file_list(fileID).name]);

    %getting lfp. rows: stimuli; columns: data (downsample to 1kHz)
    data = filtfilt(aL, bL, data');
    data=data(1:25:end);
    if length(data)>length_expected
        data=data(1:length_expected);
    elseif length(data)<length_expected
        data=[data; zeros(length_expected-length(data),1)];
    end
    lfp(:,fileID)=data;
end


save([path2save,date,'_LFP_CH',int2str(ch),'_',keyWord],'lfp');




