function [eigs negsum possum energies Matrices]=detectAxonalBreakpoint(pathToAnalysisData, patternNo)

lMovingAverage=1;
tMin = 11;
tMax = 40;
movieNos  = findMovieNos(pathToAnalysisData,patternNo);



movieNos = sort(movieNos);

dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
    movieNos(1), 99999);
firstArtifact = mean(dataTraces,1);

for movieIndex = 1:length(movieNos)

   
    dataTraces=NS_ReadPreprocessedData(pathToAnalysisData, '', 0, patternNo,...
        movieNos(movieIndex), 99999);

    subtractionMatrix = repmat(firstArtifact,[size(dataTraces,1) 1]);

    PowerMatrix=0;
    suma=0;
    for t=tMin:tMax
        meanData = mean(dataTraces(:,:,t)-subtractionMatrix(:,:,t),1);
        %meanData = mean(dataTraces(:,:,t));

       [~,EIm_view]   = ei2matrix(meanData);
       suma=suma+EIm_view;
    end
    
    suma=suma./(tMax-tMin+1);
    for t = tMin:tMax 
      
        meanData = mean(dataTraces(:,:,t)-subtractionMatrix(:,:,t),1);
        %meanData = mean(dataTraces(:,:,t));

       [~,EIm_view]   = ei2matrix(meanData);
       
       PowerMatrix=PowerMatrix+(EIm_view-suma).^2;
       
    end
    
    PowerMatrix = PowerMatrix;
    Matrices{movieIndex}= PowerMatrix;
    
    
    
A=Matrices{movieIndex};
for n=1:4;
[maxValue, linearIndexesOfMaxes] = max(A(:));
[rowsOfMaxes colsOfMaxes] = find(A == maxValue);
for i=1:length(rowsOfMaxes)
    A(rowsOfMaxes(i),colsOfMaxes(i))=NaN;
end
end
energies(movieIndex)=nansum(nansum(A));
A=A./nansum(nansum(A));
Matrices{movieIndex}=A;
x =[1:32]';
y =[1:65];
Mx= nansum(x'*A);
My = nansum(A*y');

x2=x-Mx;
y2=y-My;

suma=0;
for i=1:32
    for j=1:65
        if(~isnan(A(i,j)))
        vec=[x2(i) y2(j)]';
        mat=vec*vec';
      suma=suma+mat*A(i,j);
        end
    end
end


[a b]=eig(suma);
eigs(movieIndex,2)=max(diag(b));
eigs(movieIndex,1)=min(diag(b));

end

ratios= eigs(:,1)./(eigs(:,1)+eigs(:,2));
movingAverage=ones(1,lMovingAverage)/lMovingAverage;

ratiosFiltered=filter(movingAverage,1,ratios);
J = length(ratiosFiltered);
for i=1:J-1
    valsAfter = ratiosFiltered(i:end);
    difvalsAfter = diff(valsAfter);
    IndexNeg     = difvalsAfter<0;
    if(IndexNeg(1)==0)
            negsum(i) = 0;
    else
        if(nansum(IndexNeg)==length(IndexNeg))
            rangeDec=[i:J];
        else
            aux=find(IndexNeg==0);
            aux=aux(1);
            rangeDec=[i:i+aux(1)-1];
        end
        negsum(i)=ratiosFiltered(rangeDec(end))-ratiosFiltered(rangeDec(1));
        
    end
end


for i=1:J-1
    valsAfter = ratiosFiltered(i:end);
    difvalsAfter = diff(valsAfter);
    IndexPos     = difvalsAfter>0;
    if(IndexPos(1)==0)
            possum(i) = 0;
    else
        if(nansum(IndexPos)==length(IndexPos))
            rangeInc=[i:J];
        else
            aux=find(IndexPos==0);
            aux=aux(1);
            rangeInc=[i:i+aux(1)-1];
        end
        possum(i)=ratiosFiltered(rangeInc(end))-ratiosFiltered(rangeInc(1));
        
    end
end

end
            

