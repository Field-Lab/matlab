%Ten skrypt generuje po prostu usredniony sygnal na wszystkich 512
%elektrodach, dla wszystkich (np. 32) kolejnych movies dla danego numeru
%pattern. Numer patternu jest w nazwie. Nie ma tutaj zadnego klastrowania
%itp, po prostu usredniony sygnal ze wszystkich przebiegow.
clear
NS_GlobalConstants=NS_GenerateGlobalConstants(500);

Patterns=[1:512];
Movies=[1:2:63];

DataPath='G:\analysis\2012-09-27-4\scan_new';
EIsDataPath='D:\Home\Pawel\analysis\retina\2012-09-27-4\analysis_2013_08_06\';
Indexes=[1:50]; % in the future, this should be changed and for each pattern-movie combination, only repetitions that belong to the main cluster should be taken
EI=zeros(length(Movies),512,140);

%done: 1-162, 481-512, including all electrodes in the 4 bottom rows
for i=[163:192 449:480]
    Pattern=Patterns(i);
    EIfilename=['EI_p' num2str(Pattern)];
    fid0=fopen([EIsDataPath EIfilename],'w+');
    for j=1:length(Movies)
        [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,DataPath,0,Pattern,Movies(j),1500,Indexes);
        SDT=size(DataTraces);
        EI_m=round(reshape(mean(DataTraces),SDT(2),SDT(3)));
        EI(j,:,:)=EI_m;        
    end
    fwrite(fid0,EI,'integer*2');
    fclose(fid0);
end

break
%poni?ej - jak czyta?wynikowe pliki
fid2=fopen(EIsDataPath,'r');
a=fread(fid2,'integer*2');
b=reshape(a,32,512,TraceLength);