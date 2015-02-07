function spikeTriggeredSubUnitInput(binnedResponses,cell_resp)
 tsp = find(binnedResponses==1);
 
 nSubunits = size(cell_resp,2);
 
 figure;
 icnt=0;
 for isu=1:nSubunits
     for jsu=1:nSubunits
         if(isu>jsu)
         icnt=icnt+1;
         subplot(nSubunits*(nSubunits-1)/4,2,icnt);
         input1=cell_resp(tsp,isu);
         input2=cell_resp(tsp,jsu);
         scatter(input1,input2,0.2);
         xlabel(sprintf('SU: %d ',isu));
         ylabel(sprintf('SU: %d ',jsu));
         ax=[-5:0.1:5];
         hold on
         plot(0*ax,ax,'g');
         hold on
         plot(ax,0*ax,'g');
         end
     end
 end