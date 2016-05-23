[NumberOfMovies,NumberOfPatternsPerMovie,AllPatternsUsed,Patterns]=NS512_MoviePatterns('I:\analysis\slices\2010-09-14-0\TTX_sub\movie002',NS_GlobalConstants);

for Pattern=AllPatternsUsed(8:64)
    MoviesForPattern=find(Patterns(:,Pattern)>1);
    Movies=MoviesForPattern(2:17);    
    NS512_SpikesBasedVisualization2;    
end

break
for i=1:17
    Movie=i*8;
    NS512_SpikesBasedVisualization;
end