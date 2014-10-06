PatternNumber=33;
MovieNumber=100;

FullName_pattern=['pattern' num2str(PatternNumber) '_m' num2str(MovieNumber)];
clear Pattern;
load(FullName_pattern);
StimEl=[];
for i=2:513
    a=Pattern(i).data;
    if max(max(abs(a)))>1
        StimEl=[StimEl i-1];
    end
end
StimEl

%break;

Patterns=[1:64];
Electrode=225;
PatternsWithGivenEl=[];
for p=Patterns
    FullName_pattern=['pattern' num2str(p) '_m' num2str(MovieNumber)];
    clear Pattern;
    load(FullName_pattern);
    
    a=Pattern(Electrode+1).data;    
    if max(max(abs(a)))>1
        PatternsWithGivenEl=[PatternsWithGivenEl p];
    end
end
PatternsWithGivenEl