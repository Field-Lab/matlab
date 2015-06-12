patternsNeu1699=[477:515 556:561];

for p=1:length(patternsNeu1699)
    
[thresholdHum thresholdAlg] = fitToErfOutputAndHuman(Outputs10Jun(patternsNeu1699(p)));
thresHumNeu1699(p) = thresholdHum;
thesholdAlgNeu1699(p) = thresholdAlg;
end