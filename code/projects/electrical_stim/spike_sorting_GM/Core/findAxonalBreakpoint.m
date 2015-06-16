function [input AxonBr]  = findAxonalBreakpoint(input)
%function findAxonalBreakpoint uses Power information of different movies (defined as temporal 
%                              sums of squared traces over different electrodes) to find an axonal breakpoint
%                              defined as the last j previous the onset of axonal bundle activation
%   Input:  -input: input structure with the unnormalized input.tracesInfo.Powers
%   Output: -input: input structure with now normalized input.tracesinfo.Powers information for each condition
%                   and energies information (input.tracesInfo.energies), the normalization constants.
%           -AxonBr: The axonal breakpoint (a number between 1 and J)
% Gonzalo Mena 06/15

Powers          = input.tracesInfo.Powers;
J               = length(Powers);
numberOmit      = input.params.load.findAxon.numberOmit;
lMovingAverage  = input.params.load.findAxon.lMovingAverage;
typeEigenvalue  = input.params.load.findAxon.typeEigenvalue;
for j = 1:J
    
    A = Powers{j};
    
    for n1=1:numberOmit;
        [maxValue, linearIndexesOfMaxes] = max(A(:));
        [rowsOfMaxes colsOfMaxes] = find(A == maxValue);
        for n2=1:length(rowsOfMaxes)
            A(rowsOfMaxes(n2),colsOfMaxes(n2))=NaN;
        end
    end
    Energy = nansum(nansum(A));
    A=A./Energy;
    Powers{j} = A;
    Energies(j) = Energy;
    
    x  = [1:32]';
    y  = [1:65];
    Mx = nansum(x'*A);
    My = nansum(A*y');
    
    x2=x-Mx;
    y2=y-My;
    
    variance = 0;
    for i1=1:32
        for i2=1:65
            if(~isnan(A(i1,i2)))
                vec      = [x2(i1) y2(i2)]';
                mat      = vec*vec';
                variance = variance+mat*A(i1,i2);
            end
        end
    end
    
    
    [eigvec eigval] = eig(variance);
    eigvals(j,2)=max(diag(eigval));
    eigvals(j,1)=min(diag(eigval));
    
end

ratios          = eigvals(:,1)./(eigvals(:,1)+eigvals(:,2));
movingAverage   = ones(1,lMovingAverage)/lMovingAverage;
ratiosFiltered  = filter(movingAverage,1,ratios);

input.tracesInfo.Powers   = Powers;
input.tracesInfo.Energies = Energies;

for j = 1:J-1
    valsAfter = ratiosFiltered(j:end);
    difvalsAfter = diff(valsAfter);
    IndexNeg     = difvalsAfter<0;
    if(IndexNeg(1)==0)
        negsum(j) = 0;
    else
        if(nansum(IndexNeg)==length(IndexNeg))
            rangeDec=[j:J];
        else
            aux=find(IndexNeg==0);
            aux=aux(1);
            rangeDec=[j:j+aux(1)-1];
        end
        negsum(j)=ratiosFiltered(rangeDec(end))-ratiosFiltered(rangeDec(1));
        
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

AxonBr1 = find(negsum==min(negsum));
AxonBr2 = find(possum==max(possum));
if(typeEigenvalue==1)
    AxonBr  = AxonBr1;
else
    AxonBr  = AxonBr2;
end
if(AxonBr == J|| AxonBr ==1)
    AxonBr = [];
end
