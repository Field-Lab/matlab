function Output=DoSpikeSortingLargeScale(pathToPreparation,pathToEI,patternNo,neurons,params,varargin)


Dif1=params.global.Dif1;
DifRho2=params.global.DifRho2;
positions=params.global.positions;
x=params.array.x;


Tmax=params.global.Tmax;
nTrials=params.global.nTrials;

 

if(nargin==6)
    pathToAnalysisData=[pathToPreparation 'data00' num2str(varargin{1})];
else
    
    
     
    dirs=dir(pathToPreparation);
    cont=1;
    
    for i=1:length(dirs)
    if(length(dirs(i).name)>=4)
        aux=find(dirs(i).name(1:4)=='data');
    if(length(aux)>0)
        
        FoldersInd(cont)=str2num(dirs(i).name(5:end));
        cont=cont+1;
    
    end
    end
    end
  




[theta,rho] = cart2pol(positions(:,1)-positions(patternNo,1),positions(:,2)-positions(patternNo,2));

ind=setdiff([1:512],find(rho==0));
rho=rho(ind);


for j=1:length(listAmps)
for i=1:length(listAmps)
DifAmp1(i,j)=abs(listAmps(i)-listAmps(j));
end
end
Difs{2}=DifRho2(ind,ind)/max(max(DifRho2(ind,ind)));
Difs{3}=DifAmp1;
Difs{1}=Dif1;    
Diags{1}=[1:Tmax]'/Tmax;
Diags{2}=rho/max(rho);
Diags{3}=listAmps'/max(listAmps);
     
Difs{3}=Difs{3}/max(max(Difs{3}));
   


factp=[0 3 6 9];

for k=1:3
[Ker KerD]=evalKernels(Difs{k},Diags{k},x(factp(k)+1:factp(k+1)),[1 1 1]);
[a b]=eig(Ker);
Q{k}=a';
Qt{k}=a;
dL{k}=diag(b);
Kers{k}=Ker;
end


[TracesAll Art var0 listAmps listCurrent]=loadTracesArt(pathToAnalysisData,patternNo,Tmax,nTrials);


params.patternInfo.Q=Q;
params.patternInfo.Qt=Qt;
params.patternInfo.Kers=Kers;
params.patternInfo.dL;
params.patternInfo.ind;



[spikes Log]=SpikeSortingNoBundleNoStim(params,TracesAll)
