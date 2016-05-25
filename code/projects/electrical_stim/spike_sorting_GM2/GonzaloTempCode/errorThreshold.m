%% Error Calculations
load Results
for g=setdiff([1:7],[6])
    Outputs=RespStimBundle{g};
    neuronTemplates=templates{g};
    thresAct=0.6;
    thresHolds=[-200:-5];
    contNeurons{g}=zeros(length(Outputs{1}.neuronInfo.neuronIds),1);
    contActivation{g}=zeros(length(Outputs{1}.neuronInfo.neuronIds),1);
    contActivationT{g}=zeros(length(Outputs{1}.neuronInfo.neuronIds),length(thresHolds));
    contSameElectrode{g}=zeros(1,length(Outputs{1}.neuronInfo.neuronIds));
    indSameElectrode{g}=[];
   pat=unique(patternwr(indwr{g},1));
        
    for i=1:length(indwr{g})
        nTrialsAll{g}(i)=0;
        nTrialsAllT{g}(i,:)=zeros(1,length(thresHolds));
        neuronIds=Outputs{1}.neuronInfo.neuronIds;
        neuron=neuronIdwr(indwr{g}(i));
        
        nneu=find(neuron==neuronIds);
        contNeurons{g}(nneu)=contNeurons{g}(nneu)+1;
        patternNo=patternwr(indwr{g}(i));
        patindex=find(patternNo==pat);
        listAmps=Outputs{patindex}.stimInfo.listAmps;
        stimElec=Outputs{patindex}.stimInfo.stimElec;
        [isSame distance] = distStimSomaElectrode(neuronTemplates,stimElec);
        distances{g}(i)=distance(nneu);
        if(isSame(nneu))
            contSameElectrode{g}(nneu)=contSameElectrode{g}(nneu)+1;
            indSameElectrode{g}=[indSameElectrode{g} i];
        end
        
        pathToAnalysisData=[Outputs{patindex}.path.pathToAnalysisData ];
        name=['elecResp_n' num2str(neuron) '_p' num2str(patternNo) '.mat'];
        load([pathToAnalysisData name]);
        dif{g}{i}=NaN*zeros(length(elecResp.analysis.latencies),max(Outputs{patindex}.stimInfo.nTrials));
        nPositives{g}(i)=0;
        for j=1:length(elecResp.analysis.latencies)
            if(isempty(elecResp.analysis.latencies{j}))
                continue
            end
            nTrials=length(elecResp.analysis.latencies{j});
            dif{g}{i}(j,1:nTrials)=double((elecResp.analysis.latencies{j}>0)')-double(Outputs{patindex}.neuronInfo.spikes{nneu}(j,1:nTrials)>0);
            nPositives{g}(i)=nPositives{g}(i)+nansum(double((elecResp.analysis.latencies{j}>0)'));
            nTrialsAll{g}(i)=nTrials+nTrialsAll{g}(i);
        end
        nNegatives{g}(i)=nTrialsAll{g}(i)-nPositives{g}(i);
        error{g}(i)=nansum(nansum(abs(dif{g}{i})));
        FP{g}(i)=nansum(nansum(dif{g}{i}<0));
        FN{g}(i)=nansum(nansum(dif{g}{i}>0));
        
        contActivation{g}(nneu)=contActivation{g}(nneu)+double(max(elecResp.analysis.successRates)>thresAct);;
        current{g}(i)=max(max(abs(Outputs{patindex}.stimInfo.listAmps)));
        
        
        for p=1:length(thresHolds)
            nPositivesT{g}(i,p)=0;
            [onset onsetC]=findBundleFrompValStandAlone(Outputs{patindex}.bundle.pvals,listAmps,thresHolds(p));
            onsetsC{g}(i)=onsetC;
            if(isnan(onsetC))
                nPositivesT{g}(i,p)=nPositives{g}(i);
                nTrialsAllT{g}(i,p)=nTrialsAll{g}(i);
                errorT{g}(i,p)=error{g}(i);
                currentT{g}(i,p)=max(max(abs(Outputs{patindex}.stimInfo.listAmps)));
                FPT{g}(i,p)=FP{g}(i);
                FNT{g}(i,p)=FN{g}(i);
                
                contActivationT{g}(nneu,p)=contActivationT{g}(nneu,p)+double(max(elecResp.analysis.successRates)>thresAct);;
                 nNegativesT{g}(i,p)=nTrialsAllT{g}(i,p)-nPositivesT{g}(i,p);
            elseif(onsetC==1)
                nTrialsAllT{g}(i,p)=0;
                errorT{g}(i,p)=0;
                FPT{g}(i,p)=0;
                FNT{g}(i,p)=0;
                currentT{g}(i,p)=0;
                nPositivesT{g}(i,p)=0;
                 nNegativesT{g}(i,p)=nTrialsAllT{g}(i,p)-nPositivesT{g}(i,p);
            else
                nTrialsAllT{g}(i,p)=nansum(Outputs{patindex}.stimInfo.nTrials(1:onsetC-1));
                errorT{g}(i,p)=nansum(nansum(abs(dif{g}{i}(1:onsetC-1,:))));
                
                currentT{g}(i,p)=max(max(abs(Outputs{patindex}.stimInfo.listAmps(1:onsetC-1))));
                contActivationT{g}(nneu,p)=contActivationT{g}(nneu,p)+double(max(elecResp.analysis.successRates(1:onsetC-1))>thresAct);;
                for jj=1:onsetC-1
                    nPositivesT{g}(i,p)=nPositivesT{g}(i,p)+nansum(double((elecResp.analysis.latencies{jj}>0)'));
                    
                end
                 nNegativesT{g}(i,p)=nTrialsAllT{g}(i,p)-nPositivesT{g}(i,p);
                 
                 FPT{g}(i,p)=nansum(nansum(dif{g}{i}(1:onsetC-1,:)<0));
                 FNT{g}(i,p)=nansum(nansum(dif{g}{i}(1:onsetC-1,:)>0));
        
                 
                 
            end
            
           
        end
    end
    errorSE{g}=error{g}(indSameElectrode{g});
    errorTSE{g}=errorT{g}(indSameElectrode{g},:);
    nTrialsAllSE{g}=nTrialsAll{g}(indSameElectrode{g});
    nTrialsAllSET{g}=nTrialsAllT{g}(indSameElectrode{g},:);
end

error{9}=[];
errorT{9}=[];
FP{9}=[];
FN{9}=[];
FPT{9}=[];
FNT{9}=[];

nTrialsAll{9}=[];
nTrialsAllT{9}=[];
nPositives{9}=[];
nNegatives{9}=[];
nPositivesT{9}=[];
nNegativesT{9}=[];
contActivation{9}=[];
contActivationT{9}=[];
indwr{9}=[];
distances{9}=[];
for g=setdiff([1:7],6)
    error{9}=[error{9} error{g}];
    errorT{9}=[errorT{9};errorT{g}];
    FP{9}=[FP{9} FP{g}];
    FN{9}=[FN{9} FN{g}];
    nPositives{9}=[nPositives{9} nPositives{g}];
    nNegatives{9}=[nNegatives{9} nNegatives{g}];
    FPT{9}=[FPT{9}; FPT{g}];
    FNT{9}=[FNT{9}; FNT{g}];
    nPositivesT{9}=[nPositivesT{9};nPositivesT{g}];
    nNegativesT{9}=[nNegativesT{9};nNegativesT{g}];
    nTrialsAll{9}=[nTrialsAll{9} nTrialsAll{g}];
    nTrialsAllT{9}=[nTrialsAllT{9}; nTrialsAllT{g}];
    contActivation{9}=[contActivation{9};contActivation{g}];
    contActivationT{9}=[contActivationT{9};contActivationT{g}];
    indwr{9}=[indwr{9} indwr{g}];
    distances{9}=[distances{9} distances{g}];
end


indSameElectrode{9}=find(distances{9}==0);

%% Plot
curr=cd;
ind0=130;
for g=setdiff([1:9],[6 8])
    cd('/Users/Cybele/Google Drive/Research/spike-sorting-electrical-artifact/ResultsLargeScale/Accuracies')
    h=figure(g);
    
    set(h,'PaperUnits','inches','PaperPosition',[0 0 14 7.5])
    set(h,'position',[100   100    1500   700])
    subplot(2,2,1)
    plot(thresHolds(ind0:end),1-errorT{g}(:,ind0:end)./nTrialsAllT{g}(:,ind0:end),'linewidth',2)
    xlabel('Threshold')
    ylabel('Probability')
    title(['Individual accuracies, ' num2str(length(indwr{g})) ' Pattern-neuron pairs'])
    set(gca,'fontsize',15)
    axis([thresHolds([ind0 end]) 0.5 1])
    subplot(2,2,2)
    plot(thresHolds(:,ind0:end),nanmean(1-errorT{g}(:,ind0:end)./nTrialsAllT{g}(:,ind0:end)),'linewidth',2)
    hold on
    title('Summary measures')
    
    plot(thresHolds(:,ind0:end),nanmean(nTrialsAllT{g}(:,ind0:end)./repmat(nTrialsAll{g}',1,length(thresHolds(:,ind0:end)))),'linewidth',2)
    xlabel('Threshold')
    ylabel('Probability')
    h2=legend('Accuracy','Proportion of trials used');
    set(h2,'Position',[0.6 0.6 0.11 0.1])
    set(h2,'fontsize',10)
    set(gca,'fontsize',15)
    axis([thresHolds([ind0 end]) 0.5 1])
    
    subplot(2,2,3)
    
    % plot(tresHolds,contActivationT{g}./repmat(contActivation{g},1,length(thresHolds)),'linewidth',2)
    plot(thresHolds(:,ind0:end),contActivationT{g}(:,ind0:end)./repmat(contActivation{g},1,length(thresHolds(:,ind0:end))),'linewidth',2)
    title(['Individual neuron n=' (num2str(length(contActivation{g}))) ' activation probability'])
    xlabel('Threshold')
    ylabel('Probability')
    set(gca,'fontsize',15)
    axis([thresHolds([ind0 end]) 0 1])
    subplot(2,2,4)
    plot(thresHolds(:,ind0:end),nansum(contActivationT{g}(:,ind0:end))./nansum(contActivation{g}),'linewidth',2)
    
    xlabel('Threshold')
    ylabel('Probability')
    
    title('Summary of neuron activation')
    set(gca,'fontsize',15)
    hold on
    plot(thresHolds(:,ind0:end),nansum(contActivationT{g}(:,ind0:end)>0)./nansum(contActivation{g}>0),'linewidth',2)
    xlabel('Threshold')
    ylabel('Probability')
    h2=legend('Prob  activation neu','Prop of activated neu');
    set(h2,'Position',[0.6 0.12 0.11 0.1],'fontsize',10)
    set(gca,'fontsize',15)
    axis([thresHolds([ind0 end]) 0 1])
    print(['Preparation' num2str(g)],'-djpeg','-r100')
end
cd(curr)

for num=1:17
    SummaryPattern(Outputs{num},params,num)
end

g=1;
error1=find(1-errorT{g}(:,end)./nTrialsAllT{g}(:,end)<0.95);
%for g=1:
    
for i=1:length(indSameElectrode{g})
    pat=unique(patternwr(indwr{g},1));
    patternNo=patternwr(indwr{g}(indSameElectrode{g}(i)));
   %for patindex=26:50
    patindex=find(patternNo==pat);
    SummaryPattern(RespAllFastStim{g}{patindex},params,patindex)
    [isSame distance] = distStimSomaElectrode(templates{g},RespAllFastStim{g}{patindex}.stimInfo.stimElec);
    isSames(i,:)=isSame;
end


[a b]=sort(distances{9});

udist=unique(a);

for d=1:length(udist)
    indi=find(distances{9}==udist(d));
    errorD(d)=nansum(error{9}(indi));
    trialsD(d)=nansum(nTrialsAll{9}(indi));
end
scatter(udist,1-errorD./trialsD)
b = glmfit(udist',[errorD' trialsD'],'binomial','link','probit');


%%Error Analysis
g=2;

floor(1000*[nansum(nTrialsAll{g}) nansum(nPositives{g}) nansum(nNegatives{g}) nansum(nPositives{g})/nansum(nTrialsAll{g}) 1-nansum(error{g})./nansum(nTrialsAll{g})  (nansum(nPositives{g})-nansum(FN{g}))./nansum(nPositives{g}) (nansum(nNegatives{g})-nansum(FP{g}))./nansum(nNegatives{g})])/1000;

floor(1000*[nansum(nTrialsAllT{g}(:,189)) nansum(nPositivesT{g}(:,189)) nansum(nNegativesT{g}(:,189)) nansum(nPositivesT{g}(:,189))/nansum(nTrialsAllT{g}(:,189)) 1-nansum(errorT{g}(:,189))./nansum(nTrialsAllT{g}(:,189))  (nansum(nPositivesT{g}(:,189))-nansum(FNT{g}(:,189)))./nansum(nPositivesT{g}(:,189)) (nansum(nNegativesT{g}(:,189))-nansum(FPT{g}(:,189)))./nansum(nNegativesT{g}(:,189))])/1000



floor(1000*[nansum(nTrialsAll{g}(indSameElectrode{g})) nansum(nPositives{g}(indSameElectrode{g})) nansum(nNegatives{g}(indSameElectrode{g})) nansum(nPositives{g}(indSameElectrode{g}))/nansum(nTrialsAll{g}(indSameElectrode{g})) 1-nansum(error{g}(indSameElectrode{g}))./nansum(nTrialsAll{g}(indSameElectrode{g}))  (nansum(nPositives{g}(indSameElectrode{g}))-nansum(FN{g}(indSameElectrode{g})))./nansum(nPositives{g}(indSameElectrode{g})) (nansum(nNegatives{g}(indSameElectrode{g}))-nansum(FP{g}(indSameElectrode{g})))./nansum(nNegatives{g}(indSameElectrode{g}))])/1000;

floor(1000*[nansum(nTrialsAllT{g}( indSameElectrode{g},189)) nansum(nPositivesT{g}( indSameElectrode{g},189)) nansum(nNegativesT{g}( indSameElectrode{g},189)) nansum(nPositivesT{g}( indSameElectrode{g},189))/nansum(nTrialsAllT{g}( indSameElectrode{g},189)) 1-nansum(errorT{g}( indSameElectrode{g},189))./nansum(nTrialsAllT{g}( indSameElectrode{g},189))  (nansum(nPositivesT{g}( indSameElectrode{g},189))-nansum(FNT{g}( indSameElectrode{g},189)))./nansum(nPositivesT{g}( indSameElectrode{g},189)) (nansum(nNegativesT{g}( indSameElectrode{g},189))-nansum(FPT{g}( indSameElectrode{g},189)))./nansum(nNegativesT{g}( indSameElectrode{g},189))])/1000