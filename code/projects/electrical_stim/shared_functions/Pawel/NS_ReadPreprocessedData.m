function [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,TracesNumberLimit)
%This function shows the data on the array layout. 

if isnumeric(PatternNumber)
    PatternNumber = num2str(PatternNumber);
end

if ArtifactSubtraction==1
    FullName=[ArtifactDataPath filesep 'p' PatternNumber '_m' num2str(MovieNumber)];
    fid=fopen(FullName,'r','ieee-le');
    b=fread(fid,'int16');    
    fclose(fid);
    b0=b(1:1000);
    b1=b(1001:length(b));
    ArtifactDataTraces=reshape(b1,b0(1),b0(2),b0(3));
    %Channels=b0(3+1:3+b0(2));
    Artifact=mean(ArtifactDataTraces);
    %b0art=b0;
end

if exist([DataPath filesep 'p' PatternNumber], 'file')
    FullName=[DataPath filesep 'p' PatternNumber filesep 'p' PatternNumber '_m' num2str(MovieNumber)];
else
    FullName=[DataPath filesep 'p' PatternNumber '_m' num2str(MovieNumber)];
end

if ~exist(FullName, 'file')
    error(['Couldn''t find file: ' FullName])
end

fid=fopen(FullName,'r','ieee-le');
b=fread(fid,'int16');
fclose(fid);
b0=b(1:1000);
b1=b(1001:length(b));
DataTraces=reshape(b1,b0(1),b0(2),b0(3));
Channels=b0(3+1:3+b0(2));

if ArtifactSubtraction==1
    for i=1:b0(1)
        DataTraces(i,:,:)=DataTraces(i,:,:)-Artifact;
    end
else
    ArtifactDataTraces=DataTraces;
end 

if nargin > 5 %if TracesNumberLimit has been specified
    if TracesNumberLimit<b0(1)
        DataTraces=DataTraces(1:TracesNumberLimit,:,:);
        %ArtifactDataTraces=ArtifactDataTraces(1:TracesNumberLimit,:,:);
    end
end