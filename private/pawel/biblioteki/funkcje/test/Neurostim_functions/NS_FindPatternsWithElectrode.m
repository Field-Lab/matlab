function PatternsEl=NS_FindPatternsWithElectrode(Patterns,MovieNumber,Electrode);

PatternsEl=[];
for p=Patterns
    FullName_pattern=['pattern' num2str(p) '_m' num2str(MovieNumber)];
    clear Pattern;
    load(FullName_pattern);
    
    a=Pattern(Electrode+1).data;
    if max(max(abs(a)))>1
        PatternsEl=[PatternsEl p];
    end
end