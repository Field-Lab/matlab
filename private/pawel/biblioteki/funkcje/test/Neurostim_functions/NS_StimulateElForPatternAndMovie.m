function StimEl=NS_StimulateElForPatternAndMovie(PatternNumber,MovieNumber);

FullName_pattern=['pattern' num2str(PatternNumber) '_m' num2str(MovieNumber)];
clear Pattern;
load(FullName_pattern);
StimEl=[];
for i=2:length(Pattern)-1
    a=Pattern(i).data;
    if max(max(abs(a)))>1
        StimEl=[StimEl i-1];
    end
end