function params=InitializeArray(pathToPreparation,varargin)

load arrayPositions512

params.global.Tmax=40;
params.global.tarray=[0 [7:32]];
params.global.options=optimoptions('fminunc','Algorithm','trust-region','GradObj','on');
params.global.nTrials = 50;
params.global.x0=[2.6766    2.6729    1.5639    2.5233    1.9566  -20.7433    0.3893 32.0274];
params.global.positions=positons;

x0=params.global.x0;

tarray=params.global.tarray;
Tmax=params.global.Tmax;
nTrials=params.global.Tmax;
options=params.global.options;

path=[];
patternNo=[];


if(nargin==2)
    patternNo=varargin{1};
end

if(nargin<=2)
    
    dirs=dir(pathToAnalysisData);
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
    
elseif(nargin==3)
    FoldersInd=varargin{2};
end

for f=1:length(FoldersInd)
    pathAux=[pathToPreparation 'data00' num2str(FoldersInd(f))];
        dirs=dir(pathAux);
        
        
        
        for i=1:length(dirs)
        
        aux=find(dirs(i).name(1)=='p');
        if(length(aux)>0)
            if(nargin==1)
                if(~isequal(dirs(i).name(2),'a'))
                    patternNo=str2num(dirs(i).name(2:end));
                end
            end
            if(isequal(dirs(i).name(2:end),num2str(patternNo)))
                path=pathAux;
            end
        end
        end
end

     if(isempty(path))
         disp('Could not find useful pattern')
         params.Log{1}='Could not find useful pattern';
         return
     end
   
  pathToAnalysisData=path;   
 
Tmax=40;
    times=[1:Tmax]';
    Diags{1}=[1:Tmax]'/Tmax;
for j=1:length(times)
for i=1:length(times)

Dif1(i,j)=abs(i-j);
end
end
Dif1=Dif1/(max(max(Dif1)));
Difs{1}=Dif1;    
   
params.global.Dif1;
params.global.DifRho2;
  
for i=1:512
for j=1:512

DifRho2(i,j)=norm(positions(i,:)-positions(j,:),2);
end
end
 

[theta,rho] = cart2pol(positions(:,1)-positions(patternNo,1),positions(:,2)-positions(patternNo,2));

ind=setdiff([1:512],find(rho==0));

rho=rho(ind);


[TracesAll Art var0 listAmps listCurrent]=loadTracesArt(pathToAnalysisData,patternNo,Tmax,nTrials);%Knn2(ind,:)=squeeze(Knn(22,:,:))

for j=1:length(listAmps)
for i=1:length(listAmps)
DifAmp1(i,j)=abs(listAmps(i)-listAmps(j));
end
end
Difs{2}=DifRho2(ind,ind)/max(max(DifRho2(ind,ind)));
Difs{3}=DifAmp1;
Diags{2}=rho/max(rho);
Diags{3}=listAmps'/max(listAmps);
     
Difs{3}=Difs{3}/max(max(Difs{3}));

  f1=@(Art,x)logDetKron(Art(:,ind,:),[x(1:7) -100 -100 x(8) log(var0)],Difs,setdiff([1:10],[8 9]),[1 1 1],Diags,[3 3 3]);
 g1=@(x)f1(Art,x);

 x1 = fminunc(g1,x0,options);
 x=[x1(1:7) -100 -100 x1(end)];

 params.array.x= x;
 params.array.patternNo =patternNo;
 
 params.patternInfo.Difs=Difs;
 params.patternInfo.Diags=Diags;
