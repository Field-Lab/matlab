function CurveCorr=NS_CorrectCurve(Efficacies);
%This function takes the values in Efficacies and 
CurveCorr=Efficacies;
for i=length(Efficacies):-1:1
    if Efficacies(i)==0
        CurveCorr(i)=100;
    else
        break;
    end
end