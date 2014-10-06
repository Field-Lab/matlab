function signal=NS_SubtractArtifact_StimchipPaper(Trace,Artifact,StimulationTimings);

%STraces=size(Traces);
%SArtifactEI=size(ArtifactEI);
%LT=STraces(3)-1;
%LA=SArtifactEI(2)-1;

LT=length(Trace);
LA=length(Artifact)

for i=1:length(StimulationTimings)
    index=StimulationTimings(i)
    figure(155);
    %plot(Trace(1,index:index+LA-1));
    Trace(1,index:index+LA-1)=Trace(1,index:index+LA-1)-Artifact;
end
signal=Trace;
return;
for i=1:4000%length(StimulationTimings)
    ST0=StimulationTimings(i);   %pierwsza probka, dla ktorej nalezy wykonac odejmowanie
    Samples=[ST0:ST0+LA];        %wszystkie probki, dla ktorych....
    ST=find(Samples>0 & Samples<LT); %wszystkie z tych probek, ktore naleza do zakresu czasowego
    LST=length(ST);
    if LST>1        
        %A=reshape(ArtifactTrace(1:LST),1,LST);
        %B=reshape(SignalTrace(Samples(ST)),1,LST);        
        Traces(Samples(ST))=Traces(Samples(ST))-ArtifactEI(1:LST); % do korekty!!!
    end
end

signal=Traces;