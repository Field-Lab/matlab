function [ResSpikesE ResSpikesEJ]=substractActionPotentials(dataVecJ,ActionPotentials,E)


T=size(dataVecJ{1},2)/E;
J=size(ActionPotentials,2);
nNeurons=size(ActionPotentials,1);

for e = 1:E
ResSpikesE{e} = [];
end

for j=1:J
    I(j) = size(ActionPotentials{1,j},1);
   
    sumActionPotentials=0;
    for n = 1:nNeurons
        
        sumActionPotentials = sumActionPotentials + ActionPotentials{n,j};
    
    end
    auxResSpikes = dataVecJ{j} - sumActionPotentials;
    
    for e = 1:E
       
        ResSpikesEJ{e}{j}= auxResSpikes(:,(e-1)*T+1:e*T);
        ResSpikesE{e} =[ResSpikesE{e};reshape(ResSpikesEJ{e}{j}',T*I(j),1)];
    
    end
    

end
