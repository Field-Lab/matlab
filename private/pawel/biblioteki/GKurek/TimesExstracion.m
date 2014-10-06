clear
clc

NS_GlobalConstants = NS_GenerateGlobalConstants(512);
full_path='C:\home\Pawel\nauka\analiza\SlicesTTX\data002min009'; 

cd (full_path)
cd ..
%cd C:\home\Pawel\nauka\analiza\SlicesTTX;
%NumberOfMovies = NS512_GetNumberOfMovies(full_path);

fid0=fopen('C:\home\Pawel\nauka\analiza\SlicesTTX\movie002','r','b');
header=readMHchunk(fid0);
NumberOfMovies=header.number_of_movies;
fclose(fid0);

MoviesInfo=[];
    
for i=1:NumberOfMovies
    
    MovieData4=NS_MovieData('002',i,NS_GlobalConstants);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData4);        
    
    MoviesInfo = [MoviesInfo ; i MovieBegin];
end

paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile('C:\home\Pawel\nauka\analiza\SlicesTTX\data002minus009paramSninya\2010-09-14-0\data002min009\data002min009.params');
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\home\Pawel\nauka\analiza\SlicesTTX\data002minus009paramSninya\2010-09-14-0\data002min009\data002min009.neurons');
idList = neuronFile.getIDList(); 

NeuronID=idList(12);
spikeTimes = neuronFile.getSpikeTimes(NeuronID)'

out=[];
ignored=0;
tic
 for i=1:size(spikeTimes,2)
     
    moviesBelowTh = find(MoviesInfo(:,2) <= spikeTimes(i));
    if length(moviesBelowTh)==0
       continue
    end    
    
    MovieNo = moviesBelowTh(end);
    TimeInMovie =spikeTimes(i) - MoviesInfo(MovieNo,2);
    rest = mod(TimeInMovie,RepetPeriod);
    MovieData4=NS_MovieData('002',MovieNo,NS_GlobalConstants);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData4);
    times = MovieData(1:3:end);
    electodes = MovieData(2:3:end);
    timesBelowTh = find(times <= rest);
    if(size(timesBelowTh,1)==0) ignored =ignored+1; continue; end;
    timeNo=timesBelowTh(end);
    electrodeNo = electodes(timeNo);
    
    elem = [electrodeNo; MovieNo;  rest-times(timeNo)];
    out = [out elem];
 end
toc
el=out(1,:);
mv=out(2,:);
del=out(3,:);
unique(el)

a=hist(double(el),1000);
