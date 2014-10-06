NS_GlobalConstants=NS_GenerateGlobalConstants(500);

idList=[168 261 379 649 751 782 800 976 993 994 1026 1129 1341 1444 1787 1923 1999 2043 2164 2416]
PrimEl=[12 18 29 44 51 53 54 66 67 63 69 68 90 97 120 130 134 137 145 162]

idList=[183 261 800 858 1022 1097 1235 1264 1278 1713 2043 2222 2252 2916];
PrimEl=[18 18 54 54 66 66 83 83 40 137 149 151 195];

StimElectrodes=[127 103 174 199];

MoviePath='I:\analysis\slices\2013-12-12-3-PH\movie005';
MoviesBegins=NS512_MoviesBegins(MoviePath,NS_GlobalConstants);

RawDataPath='I:\analysis\slices\2013-12-12-3-PH\data005\';
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(RawDataPath); 

NumberOfMoviesPerElectrode=16;
MoviesForElectrode=zeros(4,16); %cztery elektrody, dwie amplitudy, osiem cz?stotliwo?ci
for i=1:length(StimElectrodes)
    indeks=1;
    StimEl=StimElectrodes(i);
    for m=1:length(MoviesBegins)-2
        MovieDataFull=NS_MovieData_GlobalPath(MoviePath,m,NS_GlobalConstants);
        StimElectrodesInTheMovie=MovieDataFull(8:3:length(MovieDataFull))
        if find(StimElectrodesInTheMovie==StimEl)
           MoviesForElectrode(i,indeks)=m;
           indeks=indeks+1;
        end
    end
end

NeuronData=zeros(16,50,10000);
for neuron=1:length(idList)                
    for i=1:length(StimElectrodes)
        StimEl=StimElectrodes(i);
        for m=1:NumberOfMoviesPerElectrode
            MovieNumber=MoviesForElectrode(i,m)
            for step=1:50                
                MoviesBegins(m)+(step-1)*10000;
                RawData=rawFile.getData(MoviesBegins(MovieNumber)+(step-1)*10000,10000)'; %the output is 65x40000 array (for 64-channel file). First index is channel number, and there is 40000 samples for each channel. The first sample is sample number 100000, as specified in the first argument.
                RepetitionData=RawData(PrimEl(neuron)+1,:);
                NeuronData(m,step,:)=RepetitionData;
            end
        end
        FullName=['D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\Matlab\RawDataPieces\ID' num2str(idList(neuron)) '_stim' num2str(StimEl)];
        fid=fopen(FullName,'wb','ieee-le');                                    
        fwrite(fid,NeuronData,'int32');
        fclose(fid);
    end
end          