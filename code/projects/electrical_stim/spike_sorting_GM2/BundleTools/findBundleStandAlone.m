function [onsets onsetsC pvals]=findBundleStandAlone(pathToPreparation,patternNo,thresHolds,Tmax,nTrials,findBundleTimes,nNeighborsExclude,varargin)
%Gonzalo Mena, 3/2016
if(nargin==4)
   FoldersNames=varargin{1};
else



    dirs=dir(pathToPreparation);
    cont=1;

    for i=1:length(dirs)
    if(length(dirs(i).name)>=4)
        aux=find(dirs(i).name(1:4)=='data');
    if(length(aux)>0)

        FoldersNames{cont}=dirs(i).name;
        cont=cont+1;

    end
    end
    end

end

for f=1:length(FoldersNames)
    pathAux=[pathToPreparation FoldersNames{f}];
        dirs=dir(pathAux);



        for i=1:length(dirs)

        aux=find(dirs(i).name(1)=='p');
        if(length(aux)>0)

            if(isequal(dirs(i).name(2:end),num2str(patternNo)))
                path=pathAux;
            end
        end
        end
end

pathToAnalysisData=path;

[TracesAll Art var0 listAmps listCurrents onset onsetC pvals]=loadTracesArt(pathToAnalysisData,patternNo,Tmax,nTrials);


[pvals]=pvalsBundleStandAlone(Art,patternNo,findBundleTimes,nNeighborsExclude);


for k=1:length(thresHolds)
    [onset onsetC]=findBundleFrompValStandAlone(pvals,listAmps,thresHolds(k));
    onsets(k)=onset;
    onsetsC(k)=onsetC;
end
