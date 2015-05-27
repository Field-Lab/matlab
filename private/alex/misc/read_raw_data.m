
samplingrate=20000;
%initiate data 
raw_data_path='/Volumes/Data/4444-44-44-4/data000';
raw_data_path='/Volumes/Data/4444-44-44-4/data001';
raw_data_path='/Volumes/Data/4444-44-44-4/data002';
raw_data_path='/Volumes/Data/4444-44-44-4/data003';

raw_data_path='/Volumes/Data/5555-55-55-5/data000';
raw_data_path='/Volumes/Data/5555-55-55-5/data001';
raw_data_path='/Volumes/Data/6666-66-66-6/data000';
raw_data_path='/Volumes/Data/7777-77-77-7/data000';
raw_data_path='/Volumes/Data/8888-88-88-8/data000';
raw_data_path='/Volumes/Data/9999-99-99-9/data000';
raw_data_path='/Volumes/Data/9999-99-99-9/data001';
raw_data_path='/Volumes/Data/2015-03-09-2/data000';
raw_data_path='/Volumes/Data/2015-03-09-2/data001';

rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(raw_data_path);

% read 1s of data starting from 0
alld = [];
for i=0:30
    i
    m1=rawFile.getData(i*samplingrate, samplingrate);
    data=m1(:,1); 
    alld = [alld; data];
end
figure
plot(alld)
tmp=diff(alld);
tmp = find(tmp<-1000);
% a = find(alld<-1000);
% a(find(diff(a)<2)+1) = [];
figure
plot(diff(tmp(10:end)))
title([raw_data_path, ', check NIL interval 100'])








for j=0:60:4400
    j
    data1=zeros(1200000,1);
    cnt=1;
    for i=j:j+59
        m1=rawFile.getData(i*samplingrate, samplingrate);
        data1(cnt:cnt+19999,1)=m1(:,1);
        cnt=cnt+20000;
    end
    
    recDuration=length(data)/samplingrate; % duration of recording in s
end

plot(data)

a = find(data<-1000)
figure
plot(diff(a))