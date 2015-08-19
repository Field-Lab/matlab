function [Output] = collapseTrialsSameCondition(data,listAmps,listStimElec,listCurrentRangesUsed)


index1=[1:size(data,1)];
cont=1;

while(~isempty(index1))
    
    amps = listAmps(index1(1),:);
    
 
    indequal = strmatch(amps,listAmps)';
        
        
        for e=1:size(data,2)
        
            dataaux{cont,e} = [];
            
            for l=indequal
            
                dataaux{cont,e} = [dataaux{cont,e}; data{l,e}];
            
            end
            
        end
   
    listAmpsAux(cont,:)                =  listAmps(index1(1),:);
    listStimElecAux(cont,:)            =  listStimElec(index1(1),:);
    listCurrentRangesUsedAux(cont,:)   =  listCurrentRangesUsed(index1(1),:);
    
    index1=setdiff(index1,indequal);
    cont=cont+1;
   
end  
Output.data                  = dataaux;
Output.listCurrentRangesUsed = listCurrentRangesUsedAux;
Output.listStimElec          = listStimElecAux;
Output.listAmps              = listAmpsAux;


         
                
                
                
             
         