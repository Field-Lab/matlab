%code to generate BARPLOTS and SCATTERPLOTS using OutputALL structure (containts 715 patterns, but 5 of them are erased, the ones
%that have a 'bad' ground truth. No weird discrepancies between Output.spikes and
%elecResp.analysis.latencies are seen, apart from one case (with lack of
%ground truth at one condition



cd('/Users/gomena/Research/EJBigData/Datasetsvisitjun15')

load OutputAll
maxi=0;

for i=1:length(OutputAll)
    recElec=OutputAll(i).tracesInfo.recElecs;
    patern=OutputAll(i).stimInfo.listStimElecs(1,:);

    crecElec=getCluster512(recElec);
    if(isempty(intersect(patern,crecElec)))
        type(i)=1;
    else
        if(isempty(intersect(patern,recElec)))
            type(i)=2;
        else
            type(i)=3;
        end
    end
end    



for i=1:length(OutputAll)
    [thresholdHum, thresholdAlg] = fitToErfOutputAndHuman(OutputAll(i));
    maxAmp(i)=max(max(abs(OutputAll(i).stimInfo.listAmps)));
    thres_humD(i)=min(maxAmp(i),max(thresholdHum,0));
    
    thres_algD(i)=min(maxAmp(i),max(thresholdAlg,0));
    if(thres_humD(i)==maxAmp(i))
        thres_humD(i)=0;
    end
    if(thres_algD(i)==maxAmp(i))
        thres_algD(i)=0;
    end
end




nNonActHum=find(thres_humD==0);
nNonActAlg=find(thres_algD==0);
nActHum=find(thres_humD>0);
nActAlg=find(thres_algD>0);
nfalsePos=length(intersect(nNonActHum,nActAlg));
nfalseNeg=length(intersect(nActHum,nNonActAlg));
ntrueNeg=length(intersect(nNonActAlg,nNonActHum));

hist(thres_humD(intersect(nActHum,nActAlg))-thres_algD(intersect(nActHum,nActAlg)))

%% Scatter plot to compare the thresholds found using the two methods
Colors(1,:)=[0 0 1];
Colors(2,:)=[255, 127, 80]/255;
Colors(3,:)=[1 0 0];
cd('/Users/gomena/Research/GIT/ChichilniskyGIT/matlab/code/projects/electrical_stim/spike_sorting_GM/Manuscript/ResultsFigures/')
titles={'50%% activation thresholds, %0.0f stim patterns\n recording electrode != stimulating electrode','50%% activation thresholds, %0.0f stim patterns\n recording electrode neighbor stimulating electrode','50%% activation thresholds, %0.0f stim patterns\n recording electrode = stimulating electrode','50%% activation thresholds'};

for i=1:3
figure; 
line(0:4,0:4,'Color','k'); hold on; 
scatter(thres_humD(type==i), thres_algD(type==i), 50, Colors(i,:),'filled');
xlabel('thresholds (\muA) using human analysis','fontsize',20); 
ylabel('thresholds (\muA) using GM algorithm','fontsize',20); 
title(titles{i},'fontsize',20)
print(['Scatter' num2str(i)],'-dpng');
print(['Scatter' num2str(i)],'-depsc2');
end
figure;
scatter(thres_humD, thres_algD, 50, Colors(type,:),'filled');
title(titles{4})
xlabel('thresholds (\muA) using human analysis','fontsize',20); 
ylabel('thresholds (\muA) using GM algorithm','fontsize',20); 
title('50% activation thresholds','fontsize',20)
print(['ScatterAll' ],'-dpng');
print(['ScatterAll' ],'-depsc2');





%% Now generate barplots

for i=1:length(OutputAll)
numTrial(i)=sum(OutputAll(i).tracesInfo.I);
end
NumTrial=sum(numTrial);
 %cases with different recording and stimulating electrodes
 
 for i=1:3
     agreement(i)=0;
     totalTrials(i)=0;
     totalHumanNegatives(i)=0;
     numFalsePositives(i) = 0;
     numTrueNegatives(i) = 0;
     totalHumanPositives(i) = 0;
     numTruePositives(i) = 0;
     numFalseNegatives(i) = 0;
     Type=find(type==i);
     for p = Type
         flag = 1;
         algorithmOutput = OutputAll(p);
         latencies = cell2mat(algorithmOutput.latencies);
         spikes = cell2mat(algorithmOutput.spikes);
         
         % Load elecResp file
         pathname = algorithmOutput.path;
         neuronId = algorithmOutput.neuronInfo.neuronIds;
         fname = ['elecResp_n' num2str(neuronId) '_p' ...
             num2str(algorithmOutput.stimInfo.patternNo) '.mat'];
         filename = fullfile(pathname,fname);
         temp = load(filename);
         elecResp = temp.elecResp;
         clear humanLat;
         humanLat = zeros(size(spikes,1), size(spikes,2));
         humanLat(:) = NaN;
         if elecResp.stimInfo.stimAmps(1,1) == elecResp.stimInfo.stimAmps(2,1)
             incSize = 2;
         else
             incSize = 1;
         end
         
         for a = 1 : size(elecResp.analysis.latencies,1)/incSize ;
             
             if elecResp.stimInfo.stimAmps(1,1) == elecResp.stimInfo.stimAmps(2,1)
                 
                 %disp('have to collapse conditions that had the same stimulation ampltiude');
                 lats1 = elecResp.analysis.latencies{2*a-1};
                 lats2 = elecResp.analysis.latencies{2*a};
                 lats=[lats1(2:end);lats2(2:end)]';
                 
                 humanLat(a,1:length(lats)) = lats;
             else
                 
                 humanLat(a,1:size(elecResp.analysis.latencies{a},1)-1) = elecResp.analysis.latencies{a}(2:end);
             end
         end
         if elecResp.stimInfo.stimAmps(1,1) == elecResp.stimInfo.stimAmps(2,1) && mod(size(elecResp.analysis.latencies,1),2)==1
             a=floor(size(elecResp.analysis.latencies,1)/incSize)+1;
             humanLat(a,1:size(elecResp.analysis.latencies{end},1)-1)=elecResp.analysis.latencies{end}(2:end);
         end
         if size(latencies,2) == (size(humanLat,2) - 1)
             humanLat(:,1) = [];
         end
         [indNanx indNany]=find(isnan(humanLat));
         humanSpikes = double(humanLat>0);
         if(~isempty(indNanx))
             for k=1:length(indNanx)
                 humanSpikes(indNanx(k),indNany(k))=-1;
             end
         end
         agreementM =  (humanSpikes == spikes);
         agreement(i) = agreement(i)+ nansum(sum(agreementM(:)));
         %     totalTrials = nansum([totalTrials  numel(humanSpikes)]);
         totalTrials(i) = totalTrials(i)+ nansum(algorithmOutput.tracesInfo.I);
         totalHumanNegatives(i) = totalHumanNegatives(i) +nansum( nansum(nansum(humanSpikes == 0)));
         numFalsePositives(i) =  numFalsePositives(i) +nansum(nansum(spikes(humanSpikes == 0)));
         numTrueNegatives(i) = numTrueNegatives(i)+nansum(sum(spikes(humanSpikes == 0) == 0));
         totalHumanPositives(i) = totalHumanPositives(i)+ nansum(nansum(nansum(humanSpikes == 1)));
         numTruePositives(i) =numTruePositives(i)+ nansum(sum(spikes(humanSpikes == 1)));
         numFalseNegatives(i) = numFalseNegatives(i)+nansum(sum(spikes(humanSpikes == 1) == 0));
     end
     
 end
 
 agreementT=sum(agreement);
 totalTrialsT=sum(totalTrials);
 numTruePositivesT=sum(numTruePositives);
 numTrueNegativesT=sum(numTrueNegatives);
 totalHumanNegativesT=sum(totalHumanNegatives);
 totalHumanPositivesT=sum(totalHumanPositives);
 
 
percentAgreement = agreement./totalTrials;
Sensitivity = numTruePositives./totalHumanPositives;
Specificity = numTrueNegatives./(totalHumanNegatives);
percentPositive= totalHumanPositives./totalTrials;

percentAgreement(4) = agreementT/totalTrialsT;
Sensitivity(4) = numTruePositivesT./totalHumanPositivesT;
Specificity(4) = numTrueNegativesT./(totalHumanNegativesT);
percentPositive(4)= totalHumanPositivesT./totalTrialsT;

all_vals = [percentAgreement;Sensitivity;...
     Specificity; percentPositive];

 
 
 for i=1:3
errors(:,i) = [1.96*sqrt(1/totalTrials(i)*percentAgreement(i)*(1-percentAgreement(i)));...
    1.96*sqrt(1/totalHumanPositives(i)*Sensitivity(i)*(1-Sensitivity(i)));...
    1.96*sqrt(1/totalHumanNegatives(i)*Specificity(i)*(1-Specificity(i)));...
    1.96*sqrt(1/totalTrials(i)*percentPositive(i)*(1-percentPositive(i)))]';
 end

 
errors(:,4) = [1.96*sqrt(1/totalTrialsT*percentAgreement(4)*(1-percentAgreement(4)));...
    1.96*sqrt(1/totalHumanPositivesT*Sensitivity(4)*(1-Sensitivity(4)));...
    1.96*sqrt(1/totalHumanNegativesT*Specificity(4)*(1-Specificity(4)));...
    1.96*sqrt(1/totalTrialsT*percentPositive(4)*(1-percentPositive(4)))]';



Colors(1,:)=[0 0 1];
Colors(2,:)=[255, 127, 80]/255;
Colors(3,:)=[1 0 0];
Colors(4,:)=[0 1 0];
%% Plots
groupnames = {' accuracy';' proportion of spikes';'sensitivity';'specificity'};
bw_colormap = Colors;
gridstatus = 'y';

figure; 
barweb([all_vals([1 4 2 3],:)], errors([1 4 2 3],:), [], ...
    groupnames, 'Spike by spike analysis of algorithm v. human', ' ',' ',bw_colormap, gridstatus);
%h=legend({'Far recording/stimulating electrodes','Close recording/stimulating electrodes', 'Same recording/stimulating electrodes','All'},'fontsize',14); 
set(gcf,'position',[  5         133        1276         571])
set(h,'position',[0.3080    0.7383    0.2257    0.1427])
set(gca,'DefaultTextFontSize',18)
set(gcf,'DefaultTextFontSize',18)
cd('/Users/gomena/Research/GIT/ChichilniskyGIT/matlab/code/projects/electrical_stim/spike_sorting_GM/Manuscript/ResultsFigures/')
print('Barplots','-dpng')
print('Barplots','-depsc2')