function SelectedMovies = NS512_GetMoviesWithPatternID(full_path,patternID)

currentWorkspace = pwd;
SelectedMovies = [];

currentWorkspace = pwd;

cd (full_path)
cd ..
if (full_path(length(full_path))=='\') shiftNTTX=1; else shiftNTTX=0; end
inputNumber = full_path(length(full_path)-2-shiftNTTX:length(full_path)-shiftNTTX);
NS_GlobalConstants = NS_GenerateGlobalConstants(512);
NumberOfMovies=NS512_GetNumberOfMovies(full_path);
for x=1:NumberOfMovies
    MovieData4=NS_MovieData(inputNumber,x,NS_GlobalConstants);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData4);
    patternIDList = MovieData(2:3:end,:);
    ind = find(patternIDList==patternID);
    
    if (size(ind,1)~=0) 
        SelectedMovies = [ SelectedMovies; x]; 
    end
end

cd (currentWorkspace);

