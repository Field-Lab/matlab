function [dataTrials]=getTrialData(raw_data_path,center_elec,time,numTrials)
%getTrialData('/Volumes/Data/2014-08-20-0/data003',center_elec,330,30)
samplingRate=20000;
%raw_data_path = '/Volumes/Data/2014-08-20-0/data003';
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(raw_data_path);
%time = 330; %seconds
display('Analysing Data');
% get Rawdata for center electrode
electrode_list=center_elec;%[0,center_elec]; 
data=zeros(time*samplingRate,length(electrode_list));
triggerData=0*data;
cnt=1;   
for j=0:1:time-1 % Depends on size of file
       
       data(cnt:cnt+1*samplingRate-1,:)=rawFile.getData(electrode_list+1,j*samplingRate, samplingRate);
       triggerData(cnt:cnt+1*samplingRate-1,:)=rawFile.getData(0+1,j*samplingRate, samplingRate);
       cnt=cnt+20000;
end

    data=data';
    data=double(data);

    
    triggerData=triggerData';
    triggerData=double(triggerData);
   

rawFile.close();
    
trialIntervalSamples=samplingRate*time/numTrials;


triggerThreshold=min(triggerData)/10;
triggerDataBinary=(triggerData<triggerThreshold);
indx=[1:samplingRate*time];
observedTriggerTimes=indx(triggerDataBinary);
correctTriggerTimes=[0:numTrials-1]*trialIntervalSamples;
trialStartIndices=[];
for iTrial=1:numTrials
[observedOffset,observedTime]=min(abs(correctTriggerTimes(iTrial)-observedTriggerTimes));
trialStartIndices=[trialStartIndices;observedTriggerTimes(observedTime)];
end



% NOTE: use only 29 samples
dataTrials=zeros(numTrials-1,trialIntervalSamples);
for iTrial=1:numTrials-1 % 30th sample less data?
    dataTrials(iTrial,:)=data(1,trialStartIndices(iTrial):trialStartIndices(iTrial)+trialIntervalSamples-1);
end

end