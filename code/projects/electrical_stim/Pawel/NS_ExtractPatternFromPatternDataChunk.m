function Pattern=NS_ExtractPatternFromPatternDataChunk(Patterns,PatternsIndexes,PatternNumber)

%extracts data for one pattern from pattern chunk data structure
%PatternsIndexes gives indeces of ends of each pattern within pattern chunk data structure

IndexEnd=PatternsIndexes(PatternNumber);
if PatternNumber==1
    IndexStart=1;
else
    IndexStart=PatternsIndexes(PatternNumber-1)+1;
end

Pattern=Patterns(IndexStart:IndexEnd);