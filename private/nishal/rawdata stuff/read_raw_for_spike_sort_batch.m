

startup

samplingRate=20000;
%initiate data
%raw_data_path='/Volumes/Data/2005-04-26-1/data006/data006000.bin';
raw_data_path = '/Volumes/Data/2012-08-09-3/data003/data003000.bin';
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(raw_data_path);

%electrode_list=[31,23,39,35,27,36,28]; % As first is trigger channel
%electrode_list =[359,355,363,367,362,354,351];
%electrode_list=[439,447,448,440,432,431,438];
%electrode_list=[495,502,503,496,487,486,494];
%electrode_list=[118,123,115,110,114,122,126];
%electrode_list=[151,142,143,152,159,158,150];
%electrode_list = [266,263,271,274,270,262,258];
 electrode_list=[1:512];
 
time = 3;%120; %seconds
data=zeros(time*samplingRate,length(electrode_list));
 
cnt=1;   
for j=0:1:time-1 % Depends on size of file
    j
       m1=rawFile.getData(j*samplingRate, samplingRate);
       data(cnt:cnt+1*samplingRate-1,:)=m1(:,electrode_list+1);
        cnt=cnt+20000;
end

rawFile.close();
    

    data=data';
    data=double(data);
   
    dt = 1/samplingRate;
    
    clear m1 raw_data_path rawFile cnt i j ans
  %  save('2005-04-26-1_nps.mat')
 %  save('../CBPSpikesortDemoPackage/example_data/2005-04-26-1_nps_7.mat');
 %save(sprintf('../CBPSpikesortDemoPackage_July9/example_data/2005-04-26_elec%d',electrode_list(1)));
 save('/Volumes/Analysis/nishal/data/nsem_1min','-v7.3')
 