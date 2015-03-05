%% Accumulate
cd('/Users/alexth/Desktop/Scripts/ONresp_scripts/')
clear

date='20130723a'
date='2013-09-03';
date='20130723b';
date='20130906';


codeWord='chirp';

mainpath=['/Users/alexth/Desktop/old_stuff/',date,'/'];
hekapath=['/Users/alexth/Desktop/old_stuff/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

path2save=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

% select heka files
clear ndFullList
file_list=[];
cnt=1;
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,codeWord, 'once'))
        file_list=[file_list i];
        
        a=regexp(heka(i).name,'ND');
        a=heka(i).name(a+2);
        if length(a)>1
            ndFullList(cnt)=str2num(a(1))+str2num(a(2));
        else
            ndFullList(cnt)=str2num(a(1));
        end
        cnt=cnt+1;
    end
end


units=dir([mainpath,'units/*.mat']);
chirp=zeros(22000,length(file_list)/3,length(units));
names=cell(length(units),1);


for cnt=1:length(units)
    load([mainpath,'units/',units(cnt).name]);
    runCnt=1;
    for i=1:3:length(file_list)        
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        conv=convolved(spikes,40,22000);
        chirp(:,runCnt,cnt)=chirp(:,runCnt,cnt)+conv(121:end-120)'; 
        runCnt=runCnt+1;
    end
    if isempty(regexp(units(cnt).name,date, 'once'))
        names{cnt}=[units(cnt).name(1),date,units(cnt).name(10:end-4)];
    else
        names{cnt}=units(cnt).name(1:end-4);
    end
end
chirp=chirp/3;

save([path2save,'chirp'],'ndFullList','chirp','names')

