function [ OutputData ] = NS512_GetData(Path, sampleChannel, StartPosition, AmountOfData)

if ( mod(AmountOfData,10000)~=0 || AmountOfData <10000 )
    error('Parameter AmountOfData in NS512_GetData function must be multiple of 10000');
end

rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(Path);
OutputData =[];    

for i=1:(AmountOfData/10000)

    RawData= rawFile.getData(StartPosition,10000)';
    OutputData = [OutputData  RawData(sampleChannel,:)];
    StartPosition = StartPosition + 10000;
end

end

