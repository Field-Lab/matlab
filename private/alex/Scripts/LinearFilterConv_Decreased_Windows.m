
clear
date='20121023';
% load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])
load(['S:\data\alexandra\MEA_data\analysis\summary_',date])
load([path2save,'protocols_',codeWord])

filter_length=500;

if length(unique(startTimes))>1
    for i=1:size(correctedProtocols,3)
        correctedProtocols(:,1,i)=correctedProtocols(:,1,i)-startTimes(i);
    end
else
    correctedProtocols(:,1,:)=correctedProtocols(:,1,:)-startTimes(1);
end

figure
for cnt=1:39
    
    load([mainpath,'units/',units(cnt).name]);
    LF=0;
    mm=1;
    m=0;
    br=0;
    for i=49:6:72
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        subplot(6,7,cnt)
        hold on
        if ~isempty(spikes)
            convSpikes=convolved(spikes,40,65000);
            convSpikes=convSpikes(121:end-120);
            br=br+basalFR;
            basalFR=mean(convSpikes(50:1950));
            
            % positiv input            
            convSpikes1=convSpikes(startTimes(i):end)-basalFR;
            convSpikes1(convSpikes1<0)=0;
%             convSpikes=-convSpikes;
            flickerTmp=flicker(:,i);
            a=zeros(length(flickerTmp)-700,700);
            for kk=1:length(flickerTmp)-700
                a(kk,700:-1:1)=flickerTmp(kk:kk+699)*convSpikes1(kk+499);
            end
            pos(mm,:)=mean(a);
%             plot(mean(a),'b','lineWidth',2)

            m=max([m max(mean(a))]);
            
            %negative input
            convSpikes1=convSpikes(startTimes(i):end)-basalFR;
            convSpikes1(convSpikes1>0)=0;
%             convSpikes1=-convSpikes1;
            flickerTmp=flicker(:,i);
            a=zeros(length(flickerTmp)-700,700);
            for kk=1:length(flickerTmp)-700
                a(kk,700:-1:1)=flickerTmp(kk:kk+699)*convSpikes1(kk+499);
            end
%             plot(mean(a),'r','lineWidth',2)
            neg(mm,:)=mean(a);
            m=max([m; max(mean(a))]);
            
            
            %common
            convSpikes1=convSpikes(startTimes(i):end)-basalFR;
%             convSpikes(convSpikes<0)=0;
%             convSpikes=-convSpikes;
            flickerTmp=flicker(:,i);
            a=zeros(length(flickerTmp)-700,700);
            for kk=1:length(flickerTmp)-700
                a(kk,700:-1:1)=flickerTmp(kk:kk+699)*convSpikes1(kk+499);
            end
%             plot(mean(a),'g')
            comm(mm,:)=mean(a);
            m=max([m max(mean(a))]);
            
            
            
            startPoint=correctedProtocols(1,1,i)+1;
            endPoint=correctedProtocols(dur,1,i);
            spikesTMP=spikes(spikes>startPoint+500&spikes<(endPoint-3000));
            spikesTMP=spikesTMP-startTimes(i);
            n=zeros(length(spikesTMP),500);
            for k=1:500
                n(:,k)=flickerTmp(spikesTMP-k+1);
            end
            LF=LF+sum(n);
            
            mm=mm+1;          
            
        end
    end
    plot(mean(pos(:,200:600)),'b','lineWidth',2)
    plot(mean(neg(:,200:600)),'r','lineWidth',2)
    plot(mean(comm(:,200:600)),'g')
    plot(LF(1:400)/max(abs(LF(1:400)))*m,'k')
    title(int2str(br/4))
    axis tight
    drawnow
end













figure
coef=1;
for cnt=1:39
    
    load([mainpath,'units/',units(cnt).name]);
    LF=0;
    mm=1;
    m=0;
    br=0;
    for i=49:2:72
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
%         subplot(6,7,cnt)
        hold on
        if ~isempty(spikes)
            convSpikes=convolved(spikes,40,65000);
            convSpikes=convSpikes(121:end-120);
            br=br+basalFR;
            basalFR=mean(convSpikes(50:1950));
            basalStd=std(convSpikes(50:1950));
            % positiv input            
            convSpikes1=convSpikes(startTimes(i):end)-basalFR;
            convSpikes1(convSpikes1<coef*basalStd)=0;
            convSpikes1(convSpikes1>coef*basalStd)=1;
            spikesTMP=find(convSpikes1(1:59000));
            spikesTMP(spikesTMP<501)=[];
            flickerTmp=flicker(:,i);
            
            n=zeros(length(spikesTMP),500);
            for k=1:500
                n(:,k)=flickerTmp(spikesTMP-k+1);
            end
            LF=sum(n);            
            plot(LF,'b')
            
            
            
            %negative input
            convSpikes1=convSpikes(startTimes(i):end)-basalFR;
            convSpikes1(convSpikes1>-coef*basalStd)=0;
            convSpikes1(convSpikes1<-coef*basalStd)=1;
            spikesTMP=find(convSpikes1(1:59000));
            spikesTMP(spikesTMP<501)=[];
            
            n=zeros(length(spikesTMP),500);
            for k=1:500
                n(:,k)=flickerTmp(spikesTMP-k+1);
            end
            LF=sum(n);
            hold on
            plot(LF,'r')
            
 
            
            startPoint=correctedProtocols(1,1,i)+1;
            endPoint=correctedProtocols(dur,1,i);
            spikesTMP=spikes(spikes>startPoint+500&spikes<(endPoint-3000));
            spikesTMP=spikesTMP-startTimes(i);
            spikesTMP(spikesTMP<501)=[];
            n=zeros(length(spikesTMP),500);
            for k=1:500
                n(:,k)=flickerTmp(spikesTMP-k+1);
            end
            LF=sum(n);
            plot(LF,'g')
            
           
            mm=mm+1;          
            
        end
    end

    title(int2str(frMean_HC(3,cnt)-frspont(3,cnt)))
    axis tight
    drawnow
end





clear
date='20120329';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])
load([path2save,'protocols_',codeWord])

filter_length=500;

if length(unique(startTimes))>1
    for i=1:size(correctedProtocols,3)
        correctedProtocols(:,1,i)=correctedProtocols(:,1,i)-startTimes(i);
    end
else
    correctedProtocols(:,1,:)=correctedProtocols(:,1,:)-startTimes(1);
end

figure
for cnt=1:59
    
    load([mainpath,'units/',units(cnt).name]);
    LF=0;
    mm=1;
    m=0;
    br=0;
    for i=9:12
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        subplot(10,6,cnt)
        hold on
        if ~isempty(spikes)
            convSpikes=convolved(spikes,40,maxLength);
            convSpikes=convSpikes(121:end-120);      
            basalFR=mean(convSpikes(50:1950));
            br=br+basalFR;
            
            % positiv input            
            convSpikes1=convSpikes(startTimes(i):end)-basalFR;
            convSpikes1(convSpikes1<0)=0;
            flickerTmp=flicker(:,i);
            a=zeros(correctedProtocols(1800,1,1)-700,700);
            for kk=1:correctedProtocols(1800,1,1)-700
                a(kk,700:-1:1)=flickerTmp(kk:kk+699)*convSpikes1(kk+499);
            end
            pos(mm,:)=mean(a);
%             plot(mean(a),'b','lineWidth',2)

            m=max([m max(mean(a))]);
            
            %negative input
            convSpikes1=convSpikes(startTimes(i):end)-basalFR;
            convSpikes1(convSpikes1>0)=0;
            flickerTmp=flicker(:,i);
            a=zeros(correctedProtocols(1800,1,1)-700,700);
            for kk=1:correctedProtocols(1800,1,1)-700
                a(kk,700:-1:1)=flickerTmp(kk:kk+699)*convSpikes1(kk+499);
            end
%             plot(mean(a),'r','lineWidth',2)
            neg(mm,:)=mean(a);
            m=max([m; max(mean(a))]);
            
            
            %common
            convSpikes1=convSpikes(startTimes(i):end)-basalFR;
            flickerTmp=flicker(:,i);
            a=zeros(correctedProtocols(1800,1,1)-700,700);
            for kk=1:correctedProtocols(1800,1,1)-700
                a(kk,700:-1:1)=flickerTmp(kk:kk+699)*convSpikes1(kk+499);
            end
%             plot(mean(a),'g')
            comm(mm,:)=mean(a);
            m=max([m max(mean(a))]);
            
            
            
            startPoint=correctedProtocols(1,1,i)+1;
            endPoint=correctedProtocols(dur,1,i);            
            spikesTMP=spikes(spikes>startPoint+500&spikes<(endPoint-3000));
            spikesTMP=spikesTMP-startTimes(i);
            spikesTMP(spikesTMP<500)=[];
            n=zeros(length(spikesTMP),500);
            for k=1:500
                n(:,k)=flickerTmp(spikesTMP-k+1);
            end
            LF=LF+sum(n);
            
            mm=mm+1;          
            
        end
    end
    plot(mean(pos(:,200:600)),'b','lineWidth',2)
    plot(mean(neg(:,200:600)),'r','lineWidth',2)
    plot(mean(comm(:,200:600)),'g')
    plot(LF(1:400)/max(abs(LF(1:400)))*m,'k')
    title(int2str(br/4))
    axis tight
    drawnow
end


           



















clear

dates=cell(15,1);
dates{1}='20130220';
dates{2}='20130220_1';
dates{3}='20130224';
dates{4}='20130225';
dates{5}='20130226';
dates{6}='20130227';
dates{7}='20130301';
dates{8}='20130301_1';
dates{9}='20130301_2';
dates{10}='20130302';
dates{11}='20130302_1';
dates{12}='20120329';
dates{13}='20120627';
dates{14}='20120714';
dates{15}='20121023';

for datesCNT=15:-1:1
    tic
    date=dates{datesCNT}    
    load(['S:\data\alexandra\MEA_data\analysis\summary_',date])
    path2save=['S:\data\alexandra\MEA_data\analysis\',date,'\'];
%     load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])
    load([path2save,'protocols_',codeWord])    
    
    if length(unique(startTimes))>1
        for i=1:size(correctedProtocols,3)
            correctedProtocols(:,1,i)=correctedProtocols(:,1,i)-startTimes(i);
        end
    else
        correctedProtocols(:,1,:)=correctedProtocols(:,1,:)-startTimes(1);
    end
    
    trialsL=min(length(file_list),192);
    posFilter=zeros(700,dim,trialsL,length(units));
    negFilter=zeros(700,dim,trialsL,length(units));
    
    for cnt=1:length(units)
        cnt
        load([mainpath,'units/',units(cnt).name]);
        
        for i=1:trialsL

            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            convSpikes=convolved(spikes,40,maxLength);
            convSpikes=convSpikes(121:end-120);
            basalFR=mean(convSpikes(50:1950));
            convSpikes=convSpikes(startTimes(i):end)-basalFR;
            flickerTmp=flicker(:,i);
            m_pos=find(convSpikes>0);
            m_neg=find(convSpikes<0);
            for t=1:dim
                if length(spikes)>10                    
                    
                    startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
                    endPoint=min(correctedProtocols(per*(t-1)+dur,1,i),length(convSpikes));                    
                    
                    
                    % positiv input
                    spikesTMP=m_pos(m_pos>(startPoint+700)&m_pos<=endPoint);
                    if length(spikesTMP)>10
                        n=zeros(length(spikesTMP),700);
                        for k=1:700
                            n(:,701-k)=flickerTmp(spikesTMP+k-500)'.*convSpikes(spikesTMP);
                        end
                        posFilter(:,t,i,cnt)=sum(n);
                    end

                    %negative input
                    spikesTMP=m_neg(m_neg>(startPoint+700)&m_neg<=endPoint);
                    if length(spikesTMP)>10
                        n=zeros(length(spikesTMP),700);
                        for k=1:700
                            n(:,701-k)=flickerTmp(spikesTMP+k-500)'.*convSpikes(spikesTMP);
                        end
                        negFilter(:,t,i,cnt)=sum(n);
                    end

                end
            end
        end
    end
    
    save([path2save,'convFilters'],'LinearFilter','posFilter','negFilter')
    toc
end








% with 1 basal std and no dim - average

clear

dates=cell(15,1);
dates{1}='20130220';
dates{2}='20130220_1';
dates{3}='20130224';
dates{4}='20130225';
dates{5}='20130226';
dates{6}='20130227';
dates{7}='20130301';
dates{8}='20130301_1';
dates{9}='20130301_2';
dates{10}='20130302';
dates{11}='20130302_1';
dates{12}='20120329';
dates{13}='20120627';
dates{14}='20120714';
dates{15}='20121023';

for datesCNT=14:-1:1
    tic
    date=dates{datesCNT}    
    load(['S:\data\alexandra\MEA_data\analysis\summary_',date])
    path2save=['S:\data\alexandra\MEA_data\analysis\',date,'\'];
%     load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])
    load([path2save,'protocols_',codeWord])     
    mainpath=['S:\data\alexandra\MEA_data\',date,'\'];
    if length(unique(startTimes))>1
        for i=1:size(correctedProtocols,3)
            correctedProtocols(:,1,i)=correctedProtocols(:,1,i)-startTimes(i);
        end
    else
        correctedProtocols(:,1,:)=correctedProtocols(:,1,:)-startTimes(1);
    end
    
    trialsL=min(length(file_list),192);
    ifDim=double(dim>1);
    posFilter_std=zeros(700,trialsL+trialsL*ifDim,length(units));
    negFilter_std=zeros(700,trialsL,length(units));
    basalFR=zeros(trialsL,length(units));
    basalFRstd=zeros(trialsL,length(units));
    metastat=zeros(trialsL+trialsL*ifDim,2); % ND, high/low contrast
    if datesCNT>5
        metastat(2:2:end,2)=1;
    end
    load([mainpath,'units/',units(1).name]);
    for i=1:size(unit{2}(file_list,1))
        tmp=cell2mat(unit{2}(file_list(i),1));
        a=regexp(tmp,'ND');
        if length(a)>1
            a=str2num(tmp(a(1)+2))+str2num(tmp(a(2)+2));
        else
            a=str2num(tmp(a+2));
        end
        metastat(i,1)=a;
    end
    
    
    
    for cnt=1:length(units)
        cnt
        load([mainpath,'units\',units(cnt).name]);
        dimCnt=1;
        for i=1:trialsL

            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            convSpikes=convolved(spikes,40,maxLength);
            convSpikes=convSpikes(121:end-120);
            basalFR(cnt,i)=mean(convSpikes(50:1950));
            basalFRstd(cnt,i)=std(convSpikes(50:1950));
            convSpikes=convSpikes(startTimes(i):end)-basalFR(cnt,i);
            flickerTmp=flicker(:,i);
            m_pos=find(convSpikes>basalFRstd(cnt,i));
            m_neg=find(convSpikes<-basalFRstd(cnt,i));
            tmp_pos=zeros(1+1*ifDim,700);
            tmp_neg=zeros(1+1*ifDim,700);
            for t=1:dim  
                cc=mod((t)-1,2)+1;
                if length(spikes)>10                    
                    
                    startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
                    endPoint=min(correctedProtocols(per*(t-1)+dur,1,i),length(convSpikes));                    
                    
                    
                    % positiv input
                    spikesTMP=m_pos(m_pos>(startPoint+700)&m_pos<=endPoint);
                    if length(spikesTMP)>10
                        n=zeros(length(spikesTMP),700);
                        for k=1:700
                            n(:,701-k)=flickerTmp(spikesTMP+k-500)'.*convSpikes(spikesTMP);
                        end
                        tmp_pos(cc,:)=tmp_pos(cc,:)+sum(n);                        
                    end

                    %negative input
                    spikesTMP=m_neg(m_neg>(startPoint+700)&m_neg<=endPoint);
                    if length(spikesTMP)>10
                        n=zeros(length(spikesTMP),700);
                        for k=1:700
                            n(:,701-k)=flickerTmp(spikesTMP+k-500)'.*convSpikes(spikesTMP);
                        end
                        tmp_neg(cc,:)=tmp_neg(cc,:)+sum(n); 

                    end

                end
            end
            posFilter_std(:,dimCnt,cnt)=tmp_pos(1,:);
            negFilter_std(:,dimCnt,cnt)=tmp_neg(1,:);            
            linFilter_noDim(:,dimCnt,cnt)=mean(LinearFilter(:,1:2:dim,i,cnt),2);
            if ifDim>0
                posFilter_std(:,dimCnt+1,cnt)=tmp_pos(2,:);
                negFilter_std(:,dimCnt+1,cnt)=tmp_neg(2,:);
                linFilter_noDim(:,dimCnt+1,cnt)=mean(LinearFilter(:,2:2:dim,i,cnt),2);
            end
            dimCnt=dimCnt+1+ifDim*1;
        end
    end
    
    save([path2save,'convFilters_std'],'posFilter_std','negFilter_std','linFilter_noDim','basalFR','basalFRstd','metastat')
%     save(['C:\daten\convFilters_std'],'posFilter_std','negFilter_std','linFilter_noDim','basalFR','basalFRstd','metastat')
   
    toc
end