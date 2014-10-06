ElectrodeNumber=512;
MovieFilePath='F:\2010-09-29-0\movie006'; %put real path here!

NS_GlobalConstants=NS_GenerateGlobalConstants(ElectrodeNumber);
[NumberOfMovies,NumberOfPatternsPerMovie,AllPatternsUsed,Patterns]=NS512_MoviePatterns(MovieFilePath,NS_GlobalConstants);
NumberOfMovies
NumberOfPatternsPerMovie
AllPatternsUsed