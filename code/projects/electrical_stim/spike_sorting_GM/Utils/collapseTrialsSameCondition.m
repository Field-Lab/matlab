function [data listAmps listStimElec  listCurrentRangesUsed] = collapseTrialsSameCondition(data,listAmps,listStimElec,listCurrentRangesUsed)


index=[1:size(data,1)];
cont=1;

while(~isempty(index))
    
    amps = listAmps(index(1),:);
    
 
    indequal = strmatch(amps,listAmps)';
        
        
        for e=1:size(data,2)
        
            dataaux{cont,e} = [];
            
            for l=indequal
            
                dataaux{cont,e} = [dataaux{cont,e}; data{l,e}];
            
            end
            
        end
   
    listAmpsAux(cont,:)                =  listAmps(index(1),:);
    listStimElecAux(cont,:)            =  listStimElec(index(1),:);
    listCurrentRangesUsedAux(cont,:)   =  listCurrentRangesUsed(index(1),:);
    
    index=setdiff(index,indequal);
    cont=cont+1;
   
end  
data                  = dataaux;
listCurrentRangesUsed = listCurrentRangesUsedAux;
listStimElec          = listStimElecAux;
listAmps              = listAmpsAux;


         
                
                
                
             
         