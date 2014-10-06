function Traces=NS512_SubtractArtifact(Traces,Artifact);

SDT=size(Traces);
for i=1:SDT(1)
    Traces(i,:,:)=Traces(i,:,:)-Artifact;
end