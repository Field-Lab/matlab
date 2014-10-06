function NewTraces=NS_PredictAndSubtractArtifact(Traces,Amplitude,Art1,Amp1,Art2,Amp2);

NewTraces=Traces;
Artifact=NS_PredictArtifact(Amplitude,Art1,Amp1,Art2,Amp2);
size(Artifact)
SNT=size(Traces);
for i=1:SNT(1)
    NewTraces(i,:,:)=Traces(i,:,:)-Artifact;
end