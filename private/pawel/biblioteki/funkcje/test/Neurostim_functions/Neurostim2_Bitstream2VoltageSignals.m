function[clk,data]=Neurostim2_Bitstream2VoltageSignals(Bitstream,SamplingPeriod,FilesPath,TriggerRecordDelay,TriggerConnectDelay);

clk=zeros(2*length(Bitstream),2);
max([1:length(clk)]*SamplingPeriod/2)
clk(:,1)=[1:length(clk)]*SamplingPeriod/2;
clk(2:2:length(clk),2)=3.3;

data=zeros(length(Bitstream),2);
data(:,1)=[1:length(Bitstream)]*SamplingPeriod;
max([1:length(Bitstream)]*SamplingPeriod)
data(find(Bitstream==1),2)=3.3;

figure(5)
plot(clk(:,1),clk(:,2),'bd-',data(:,1),data(:,2),'gd-');

fid=fopen([FilesPath '\clk.txt'],'w');
fprintf(fid,'%e  %e\n',clk');
fclose(fid);

fid=fopen([FilesPath '\data.txt'],'w');
fprintf(fid,'%e  %e\n',data');
fclose(fid);