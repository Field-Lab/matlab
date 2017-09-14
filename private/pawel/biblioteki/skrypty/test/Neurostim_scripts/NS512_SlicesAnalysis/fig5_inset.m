GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\NeuronAnalysis\2010-09-14-0\data002\GaussesFiles\';

NeuronFilePath = 'D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002\data002.neurons';
%MovieFilePath='C:\home\Pawel\nauka\512stim_paper\MovieFiles\2010-09-14-0\movie002';
DuplicatesFilePath='D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002_duplicates.txt';
PrimaryNeurons=NS512_IdentifyPrimaryNeurons_v2(NeuronFilePath,DuplicatesFilePath);

fid=fopen([GaussesFilesPath 'StrangeStim'],'r');
a=fread(fid,'double');
StrangeStim=reshape(a,length(a)/4,4);
fclose(fid);

s=0;
IloscSpikow=[];
figure(1)
clf
for n=1:length(PrimaryNeurons)
    NeuronID=PrimaryNeurons(n);
    fid=fopen([GaussesFilesPath 'AllGausses_n' num2str(PrimaryNeurons(n))],'r');
    a=fread(fid,'double');
    length(a);
    fclose(fid);
    
    b0=reshape(a,length(a)/6,6);  
    GoodFit=find(b0(:,6)<=0.1);
    b=b0(GoodFit,:);        
    
    sb=size(b);
    if sb(1)>0
        %plot(b(:,4),b(:,5),'bd');
        hold on
        b1=b;
        DobrePatterny=b(:,1);
        DobrePatternyUnique=unique(DobrePatterny);
        DaneDlaNeuronu=find(StrangeStim(:,1)==NeuronID)';
        for p=1:length(DobrePatternyUnique)
            Dane=find(StrangeStim(DaneDlaNeuronu,2)==DobrePatternyUnique(p));
            IloscSpikow=[IloscSpikow' StrangeStim(DaneDlaNeuronu(Dane),[3:4])']';            
            %tau=
        end
        
        %Patterny=StrangeStim(DaneDlaNeuronu,2)'
        
        %GaussyDlaNeuronu=find(b0(:,2)==NeuronID)
    end                    
    s=s+sb(1);  
end
figure(2)
x=[1:length(IloscSpikow)]
plot(x,IloscSpikow(:,1)./IloscSpikow(:,2))
grid on