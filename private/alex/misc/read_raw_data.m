
%initiate data 
raw_data_path='/Volumes/Data/6666-66-66-6/data000';
raw_data_path='/Volumes/Data/2015-03-09-2/data000';

raw_data_path='/Volumes/Data/8888-88-88-8/data002';

rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(raw_data_path);
samplingrate=20000;
% read 1s of data starting from 0
all_diffs = [];
part_d = []
for i=0:100:7200
    i
    for j=i:i+99        
        m1=rawFile.getData(j*samplingrate, samplingrate);        
        %     data=find(diff(m1(:,1))<-1000)+i*20000;
        %     alld = [alld; data];
        part_d = [part_d; m1(:,1)];
    end
    tmp = diff(double(part_d));
    tmp = find(tmp<-1000);
    tmp = tmp/20; % in ms
    tmp = diff((tmp));
    all_diffs = [all_diffs; tmp];
    part_d = [];
end

figure
plot(round(all_diffs(1:end)))
plot(new_time_axis(1:end-22),all_diffs)


load('/Users/alexth/Desktop/long_test_minimized.mat')
a = time_stamps{end};
new_time_axis = linspace(0,a(end), length(a));
diff_ts = diff(a);
figure
plot(new_time_axis(1:end-1), diff_ts)
xlabel('time, s')

figure
plot(diff_ts(1000:30000))
b = all_diffs(1:70000);
c = diff_ts(21:70020);

d = corr(c',b);

figure
plot(c*1000)
hold on
plot(b)


time_axis = linspace(0,tmp(end)/1000, length(tmp));
figure
plot(time_axis(1:end-1),diff(tmp));
xlabel('time, s')

figure
a = diff(tmp);
hist(a(a<11))



load('/Users/alexth/Desktop/saved_stim/04-Aug-2015/data001/time_stamps/stimulus_1.mat')
new_time_axis = linspace(0,time_stamps(end), length(time_stamps));
diff_ts = diff(time_stamps);
figure
plot(new_time_axis(1:end-1), diff_ts)
xlabel('time, s')

figure
plot(new_time_axis(new_time_axis<time_axis(end)),diff_ts(new_time_axis<time_axis(end)));
xlabel('time, s')


sum(diff_ts>0.010)





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