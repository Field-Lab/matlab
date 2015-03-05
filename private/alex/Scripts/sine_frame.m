clear
cd('/mnt/muench_data/user/alexandra/scripts')
date='20130301';
codeWord='sine';


mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

if ~exist([mainpath 'LinearFilters/'],'dir')    
    mkdir([mainpath, 'LinearFilters']);
end
path2save=[mainpath 'LinearFilters/'];

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
abused_accum=round(abused_accum);


flicker=zeros(size(abused_accum,1)+100,length(file_list));
startingPoints=zeros(1,length(file_list));
maxLength=zeros(1,length(file_list));

for i=1:length(file_list)
    flips=diff(abused_accum(1:end-2,i));
    startingPoints(i)=abused_accum(1,i);
    flips(flips>30)=2;
    flips(flips>2)=1;
    pxls=(protocol_accum(2:end-2,2,i)-30)/30;
    
    k=find(flips==2);
    cnt=0;
    for j=k
        pxls=[pxls(1:j+cnt);pxls(j+cnt);pxls(j+1+cnt:end)];
        cnt=cnt+1;
    end
    
    flicker(1:length(pxls),i)=pxls;
    maxLength(i)=length(pxls);
end
flicker=flicker(1:min(maxLength),:);

save([path2save,date,'_flicker_',codeWord,'_frame'],'flicker','startingPoints','maxLength')
