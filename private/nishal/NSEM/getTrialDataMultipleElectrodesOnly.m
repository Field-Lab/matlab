function [allElecData]=getTrialDataMultipleElectrodes(raw_data_path,center_elec,time,elecList)

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

end