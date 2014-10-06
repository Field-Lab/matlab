function Pattern=NS_ExtractPatternFromPatternDataChunk(Patterns,PatternsIndexes,PatternNumber);

IndexEnd=PatternsIndexes(PatternNumber);
if PatternNumber==1
    IndexStart=1;
else
    IndexStart=PatternsIndexes(PatternNumber-1)+1;
end

Pattern=Patterns(IndexStart:IndexEnd);