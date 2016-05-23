NeuronFilePath = 'D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002\data002.neurons';
DuplicatesFilePath='D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002_duplicates.txt';
PrimaryNeurons=NS512_IdentifyPrimaryNeurons_v2(NeuronFilePath,DuplicatesFilePath);

figure(1)
%ciekawe przypadki: 232, 
for n=201%1:length(PrimaryNeurons)
    fid=fopen(['D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002_Matlab\SpikeTimesCombined\ID=' num2str(PrimaryNeurons(n)) 'c'],'r','ieee-le'); 
    a=fread(fid,'int32');
    fclose(fid);
    l=length(a);
    b0=reshape(a,5,l/5);
    pst=b0(4,:);
    b=b0(:,find(pst>20));
    %b=b0;
        
    p1=b(3,:);
    PatternsStim=unique(p1);
    for i=2:length(unique(p1))
        ns2(i)=length(find(b(3,:)==PatternsStim(i)));
        subplot(4,4,i)
        hist(b(4,find(b(3,:)==PatternsStim(i))))
    end
    length(b)/64
    plot(ns2)
    ile(n)=length(find(ns2>=50));
end