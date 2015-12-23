function [dataTrials,allElecData,trialStartIndices]=getTrialDataMultipleElectrodes(raw_data_path,center_elec,time,numTrials,elecList)
samplingRate=20000;
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(raw_data_path);
display('Analysing Data');


% get Rawdata for center electrode

data=zeros(time*samplingRate,1);
triggerData=0*data;
allElecData=zeros(time*samplingRate,length(elecList));
cnt=1;   
for j=0:1:time-1 % Depends on size of file
        j
        m1=rawFile.getData(j*samplingRate, samplingRate);

       data(cnt:cnt+1*samplingRate-1,:)=   m1(:,center_elec+1); 
       triggerData(cnt:cnt+1*samplingRate-1,:)=m1(:,0+1);
       allElecData(cnt:cnt+1*samplingRate-1,:)=m1(:,elecList+1); 
       cnt=cnt+20000;
end

    data=data';
    data=double(data);
    
    triggerData=triggerData';
    triggerData=double(triggerData);
    
    allElecData=allElecData';
    allElecData=double(allElecData);

rawFile.close();
    
trialIntervalSamples=samplingRate*time/numTrials;


triggerThreshold=min(triggerData)/10;
triggerDataBinary=(triggerData<triggerThreshold);
indx=[1:samplingRate*time];
observedTriggerTimes=indx(triggerDataBinary);
correctTriggerTimes=[0:numTrials-1]*trialIntervalSamples;
trialStartIndices=[];
for iTrial=1:numTrials
[~,observedTime]=min(abs(correctTriggerTimes(iTrial)-observedTriggerTimes));
trialStartIndices=[trialStartIndices;observedTriggerTimes(observedTime)];
end



% NOTE: use only 29 samples
dataTrials=zeros(numTrials-1,trialIntervalSamples);
for iTrial=1:numTrials-1 % 30th sample less data?
    dataTrials(iTrial,:)=data(1,trialStartIndices(iTrial):trialStartIndices(iTrial)+trialIntervalSamples-1);
end

end