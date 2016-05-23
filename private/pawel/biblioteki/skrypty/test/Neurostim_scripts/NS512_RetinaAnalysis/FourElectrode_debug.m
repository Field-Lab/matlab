%clear

NS_GlobalConstants=NS_GenerateGlobalConstants(500);
cd D:\Home\Data\retina\LaurenGrosberg;
MovieNumber=2;
ElectrodeNumber=11;
PatternNumber=445;

MovieData=NS_MovieData('005',2,NS_GlobalConstants);
PatternsUsed=MovieData(8:3:length(MovieData));
if find(PatternsUsed==PatternNumber)
    find(PatternsUsed==PatternNumber);
else
    warning('this patterns is not used in this movie')
end
    
%[patterns_out,PatternsIndexes,Status]=ReadPatternDataChunk('D:\Home\Data\retina\LaurenGrosberg\pattern005',2,NS_GlobalConstants);

PatternIndex=PatternsIndexes(PatternNumber)
PatternIndexPrevious=PatternsIndexes(PatternNumber-1)+1

for i=PatternIndexPrevious:PatternIndex
    channel=patterns_out(i).channel
    data=patterns_out(i).data
end



%break
for i=1:16239
    c1(i)=patterns_out(i).channel;
end

plot(c1,'bd-')

p1=find(c1==ElectrodeNumber);
