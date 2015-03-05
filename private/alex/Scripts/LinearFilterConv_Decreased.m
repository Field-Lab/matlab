cd('/mnt/muench_data/user/alexandra/scripts')

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

for datesCNT=13:15
    tic
    date=dates{datesCNT}    
    
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])
    load([path2save,'protocols_',codeWord])    
    
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
    tt=1;
    for i=1:size(unit{2}(file_list,1),1)
        tmp=cell2mat(unit{2}(file_list(i),1));
        a=regexp(tmp,'ND');
        if length(a)>1
            a=str2num(tmp(a(1)+2))+str2num(tmp(a(2)+2));
        else
            a=str2num(tmp(a+2));
        end
        if ifDim>0
            metastat(tt,1)=a;
            metastat(tt+1,1)=a;
            tt=tt+2;
        else
            metastat(i,1)=a;
        end
    end
    
    
    
    for cnt=1:length(units)
        cnt
        load([mainpath,'units/',units(cnt).name]);
        dimCnt=1;
        for i=1:trialsL

            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            convSpikes=convolved(spikes,40,maxLength);
            convSpikes=convSpikes(121:end-120);
            basalFR(i,cnt)=mean(convSpikes(50:1950));
            basalFRstd(i,cnt)=std(convSpikes(50:1950));
            convSpikes=convSpikes(startTimes(i):end)-basalFR(i,cnt);
            flickerTmp=flicker(:,i);
            m_pos=find(convSpikes>basalFRstd(i,cnt));
            m_neg=find(convSpikes<-basalFRstd(i,cnt));
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
    
    save([path2save,'convFilters_1std'],'posFilter_std','negFilter_std','linFilter_noDim','basalFR','basalFRstd','metastat')
    toc
end





cd('/mnt/muench_data/user/alexandra/scripts')

% estimate firing rate parameters

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

for datesCNT=1:15
    tic
    date=dates{datesCNT}    
    
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])
    load([path2save,'protocols_',codeWord])    
    
    if length(unique(startTimes))>1
        for i=1:size(correctedProtocols,3)
            correctedProtocols(:,1,i)=correctedProtocols(:,1,i)-startTimes(i);
        end
    else
        correctedProtocols(:,1,:)=correctedProtocols(:,1,:)-startTimes(1);
    end
    
    trialsL=min(length(file_list),192);
    ifDim=double(dim>1);
    pos_mean=zeros(trialsL+trialsL*ifDim,length(units));
    neg_mean=pos_mean;
    pos_std=pos_mean;
    neg_std=pos_mean;
    basalFR=zeros(trialsL+trialsL*ifDim,length(units));
    basalFRstd=zeros(trialsL+trialsL*ifDim,length(units));
    metastat=zeros(trialsL+trialsL*ifDim,2); % ND, high/low contrast
    if datesCNT>5
        metastat(2:2:end,2)=1;
    end
    load([mainpath,'units/',units(1).name]);
    tt=1;
    for i=1:size(unit{2}(file_list,1))
        tmp=cell2mat(unit{2}(file_list(i),1));
        a=regexp(tmp,'ND');
        if length(a)>1
            a=str2num(tmp(a(1)+2))+str2num(tmp(a(2)+2));
        else
            a=str2num(tmp(a+2));
        end
        if ifDim>0
            metastat(tt,1)=a;
            metastat(tt+1,1)=a;
            tt=tt+2;
        else
            metastat(i,1)=a;
        end
    end
    
    
    
    for cnt=1:length(units)
        cnt
        load([mainpath,'units/',units(cnt).name]);
        dimCnt=1;
        for i=1:trialsL

            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            convSpikes=convolved(spikes,40,maxLength);
            convSpikes=convSpikes(121:end-120);
            basalFR(dimCnt,cnt)=mean(convSpikes(50:1950));
            basalFRstd(dimCnt,cnt)=std(convSpikes(50:1950));
            convSpikes=convSpikes(startTimes(i):end)-basalFR(i,cnt);

            m_pos=find(convSpikes>2*basalFRstd(i,cnt));
            m_neg=find(convSpikes<-2*basalFRstd(i,cnt));
            meanPos=zeros(1+1*ifDim,1);
            meanNeg=zeros(1+1*ifDim,1);
            stdPos=zeros(1+1*ifDim,1);
            stdNeg=zeros(1+1*ifDim,1);
            for t=1:dim  
                cc=mod((t)-1,2)+1;
                if length(spikes)>10                    
                    
                    startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
                    endPoint=min(correctedProtocols(per*(t-1)+dur,1,i),length(convSpikes));                    
                    
                    
                    % positiv input
                    meanPos(cc)=meanPos(cc)+nanmean(convSpikes(m_pos(m_pos>(startPoint+700)&m_pos<=endPoint)));
                    stdPos(cc)=stdPos(cc)+nanstd(convSpikes(m_pos(m_pos>(startPoint+700)&m_pos<=endPoint)));

                    %negative input
                    meanNeg(cc)=meanNeg(cc)+nanmean(convSpikes(m_neg(m_neg>(startPoint+700)&m_neg<=endPoint)));
                    stdNeg(cc)=stdNeg(cc)+nanstd(convSpikes(m_neg(m_neg>(startPoint+700)&m_neg<=endPoint)));

                end
            end
            pos_mean(dimCnt,cnt)=meanPos(1)/ceil(dim/2);
            pos_std(dimCnt,cnt)=stdPos(1)/ceil(dim/2);
            neg_mean(dimCnt,cnt)=meanNeg(1)/ceil(dim/2);
            neg_std(dimCnt,cnt)=stdNeg(1)/ceil(dim/2);            
            if ifDim>0
                pos_mean(dimCnt+1,cnt)=meanPos(2)/ceil(dim/2);
                pos_std(dimCnt+1,cnt)=stdPos(2)/ceil(dim/2);
                neg_mean(dimCnt+1,cnt)=meanNeg(2)/ceil(dim/2);
                neg_std(dimCnt+1,cnt)=stdNeg(2)/ceil(dim/2);
                basalFR(dimCnt+1,cnt)=basalFR(dimCnt,cnt);
                basalFRstd(dimCnt+1,cnt)=basalFRstd(dimCnt,cnt);
            end
            dimCnt=dimCnt+1+ifDim*1;
        end
    end
    
    save([path2save,'firingRate_pos_neg'],'pos_mean','pos_std','neg_mean','neg_std','basalFR','basalFRstd','metastat')
    toc
end





cd('/mnt/muench_data/user/alexandra/scripts')

% estimate firing rate parameters  -  numbers

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

for datesCNT=1:15
    tic
    date=dates{datesCNT}    
    
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])
    load([path2save,'protocols_',codeWord])    
    
    if length(unique(startTimes))>1
        for i=1:size(correctedProtocols,3)
            correctedProtocols(:,1,i)=correctedProtocols(:,1,i)-startTimes(i);
        end
    else
        correctedProtocols(:,1,:)=correctedProtocols(:,1,:)-startTimes(1);
    end
    
    trialsL=min(length(file_list),192);
    ifDim=double(dim>1);
    pos_mean=zeros(trialsL+trialsL*ifDim,length(units));
    neg_mean=pos_mean;
    pos_std=pos_mean;
    neg_std=pos_mean;
    basalFR=zeros(trialsL+trialsL*ifDim,length(units));
    basalFRstd=zeros(trialsL+trialsL*ifDim,length(units));
    metastat=zeros(trialsL+trialsL*ifDim,2); % ND, high/low contrast
    if datesCNT>5
        metastat(2:2:end,2)=1;
    end
    load([mainpath,'units/',units(1).name]);
    tt=1;
    for i=1:size(unit{2}(file_list,1))
        tmp=cell2mat(unit{2}(file_list(i),1));
        a=regexp(tmp,'ND');
        if length(a)>1
            a=str2num(tmp(a(1)+2))+str2num(tmp(a(2)+2));
        else
            a=str2num(tmp(a+2));
        end
        if ifDim>0
            metastat(tt,1)=a;
            metastat(tt+1,1)=a;
            tt=tt+2;
        else
            metastat(i,1)=a;
        end
    end
    
    
    
    for cnt=1:length(units)
        cnt
        load([mainpath,'units/',units(cnt).name]);
        dimCnt=1;
        for i=1:trialsL

            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            convSpikes=convolved(spikes,40,maxLength);
            convSpikes=convSpikes(121:end-120);
            basalFR(dimCnt,cnt)=mean(convSpikes(50:1950));
            basalFRstd(dimCnt,cnt)=std(convSpikes(50:1950));
            convSpikes=convSpikes(startTimes(i):end)-basalFR(i,cnt);

            m_pos=find(convSpikes>2*basalFRstd(i,cnt));
            m_neg=find(convSpikes<-2*basalFRstd(i,cnt));
            meanPos=zeros(1+1*ifDim,1);
            meanNeg=zeros(1+1*ifDim,1);
            stdPos=zeros(1+1*ifDim,1);
            stdNeg=zeros(1+1*ifDim,1);
            for t=1:dim  
                cc=mod((t)-1,2)+1;
                if length(spikes)>10                    
                    
                    startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
                    endPoint=min(correctedProtocols(per*(t-1)+dur,1,i),length(convSpikes));                    
                    
                    
                    % positiv input
                    meanPos(cc)=meanPos(cc)+sum(m_pos>(startPoint+700)&m_pos<=endPoint);
                    stdPos(cc)=stdPos(cc)+length((startPoint+700):endPoint);

                    %negative input
                    meanNeg(cc)=meanNeg(cc)+sum(m_neg>(startPoint+700)&m_neg<=endPoint);
                    stdNeg(cc)=stdNeg(cc)+length((startPoint+700):endPoint);

                end
            end
            pos_mean(dimCnt,cnt)=meanPos(1)/ceil(dim/2);
            pos_std(dimCnt,cnt)=stdPos(1)/ceil(dim/2);
            neg_mean(dimCnt,cnt)=meanNeg(1)/ceil(dim/2);
            neg_std(dimCnt,cnt)=stdNeg(1)/ceil(dim/2);            
            if ifDim>0
                pos_mean(dimCnt+1,cnt)=meanPos(2)/ceil(dim/2);
                pos_std(dimCnt+1,cnt)=stdPos(2)/ceil(dim/2);
                neg_mean(dimCnt+1,cnt)=meanNeg(2)/ceil(dim/2);
                neg_std(dimCnt+1,cnt)=stdNeg(2)/ceil(dim/2);
                basalFR(dimCnt+1,cnt)=basalFR(dimCnt,cnt);
                basalFRstd(dimCnt+1,cnt)=basalFRstd(dimCnt,cnt);
            end
            dimCnt=dimCnt+1+ifDim*1;
        end
    end
    
    save([path2save,'firingRate_pos_neg_number'],'pos_mean','pos_std','neg_mean','neg_std','basalFR','basalFRstd','metastat')
    toc
end


date='20121023'

load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date],'path2save')
load([path2save,'firingRate_pos_neg'])

figure
for cnt=1:39
    subplot(6,7,cnt)
    hold on
    cc=1;
    for i=1:24:24*8
        a(cc,1)=nanmean(pos_mean(i:i+23,cnt));
        
        a(cc,2)=nanmean(neg_mean(i:i+23,cnt));
        
        a(cc,3)=nanmean(basalFR(i:i+23,cnt));
        cc=cc+1;
    end
    plot(a)
    plot(3,a(3,:),'*')
    line([0 9],[0 0],'color','k')
end




date='20121023'

load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date],'path2save')
load([path2save,'firingRate_pos_neg_number'])

figure
for cnt=1:39
    subplot(6,7,cnt)
    hold on
    cc=1;
    for i=1:24:24*8
        a(cc,1)=nanmean(pos_mean(i:i+23,cnt))/nanmean(pos_std(i:i+23,cnt));
        
        a(cc,3)=nanmean(neg_mean(i:i+23,cnt))/nanmean(neg_std(i:i+23,cnt));
        
        a(cc,2)=nanmean(basalFR(i:i+23,cnt))/1901;
        cc=cc+1;
    end
    plot(a)
    plot(3,a(3,:),'*')
    line([0 9],[0 0],'color','k')
    
end














cd('/mnt/muench_data/user/alexandra/scripts')

% estimate firing rate parameters  -  numbers

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

for datesCNT=15
    tic
    date=dates{datesCNT}    
    
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])
    load([path2save,'protocols_',codeWord])    
    
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
    tt=1;
    for i=1:size(unit{2}(file_list,1))
        tmp=cell2mat(unit{2}(file_list(i),1));
        a=regexp(tmp,'ND');
        if length(a)>1
            a=str2num(tmp(a(1)+2))+str2num(tmp(a(2)+2));
        else
            a=str2num(tmp(a+2));
        end
        if ifDim>0
            metastat(tt,1)=a;
            metastat(tt+1,1)=a;
            tt=tt+2;
        else
            metastat(i,1)=a;
        end
    end
    
    
    
    for cnt=1:length(units)
        cnt=6
        load([mainpath,'units/',units(cnt).name]);
        dimCnt=1;
        
        
        filters=zeros(500,10);
        count=zeros(1,10);
        for i=49:2:72%:trialsL
            
            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            spikes=convolved(spikes,40,62000);
            spikes=spikes(121:end-120);
%             bfr=spikes(50:1950);
            spikes(1:startTimes(i))=[];
            spikes=diff(spikes);
%             
%             if bfr<10
%                 bfr=mean(spikes);
%             end
            bfr=spikes;
            
%             m=bfr/5;
%             m=0:m:(10*m);
%             m(end)=1000;
            m=min(bfr):(max(bfr)-min(bfr))/10:max(bfr);

            flickerTmp=flicker(:,i);
%             filters=zeros(350,15);

            for t=1:dim  
                cc=mod((t)-1,2)+1;
                if length(spikes)>10                    
                    
                    startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
                    endPoint=correctedProtocols(end,1,i)+1;
                    
                    for k=startPoint:5:endPoint-1000
                        a=mean(spikes(350+k:400+k));
                        
                        filters(:,find(m>a,1)-1)=filters(:,find(m>a,1)-1)+flickerTmp(k:(499+k));
                        count(find(m>a,1)-1)=count(find(m>a,1)-1)+1;
                    end

                    
                    
%                     for k=startPoint:35:endPoint-350
%                         a=sum(spikes>=(350+k)&spikes<(385+k));
%                         filters(:,a+1)=filters(:,a+1)+flickerTmp(k:(349+k));
%                         count(a+1)=count(a+1)+1;
%                     end
                             
                end
            end
        end
        figure
        for j=1:10
            subplot(5,2,j)
            hold on
            plot(filters(:,j)/count(j))
            title(['ND6, ', int2str(cnt)])
        end
        
        filters=zeros(500,10);
        count=zeros(1,10);
        for i=50:2:72%:trialsL
            
            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            spikes=convolved(spikes,40,62000);
            %             spikes=spikes-startTimes(i);
            spikes=spikes(121:end-120);
            bfr=mean(spikes(50:1950));
            spikes(1:startTimes(i))=[];
            
            if bfr<10
                bfr=mean(spikes);
            end
            m=bfr/5;
            m=0:m:(10*m);
            m(end)=1000;
            %             m=min(spikes):m:max(spikes);

            flickerTmp=flicker(:,i);
%             filters=zeros(350,15);

            for t=1:dim  
                cc=mod((t)-1,2)+1;
                if length(spikes)>10                    
                    
                    startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
                    endPoint=correctedProtocols(end,1,i)+1;
                    
                    for k=startPoint:5:endPoint-1000
                        a=mean(spikes(350+k:400+k));
                        
                        filters(:,find(m>a,1)-1)=filters(:,find(m>a,1)-1)+flickerTmp(k:(499+k));
                        count(find(m>a,1)-1)=count(find(m>a,1)-1)+1;
                    end

                    
                    
%                     for k=startPoint:35:endPoint-350
%                         a=sum(spikes>=(350+k)&spikes<(385+k));
%                         filters(:,a+1)=filters(:,a+1)+flickerTmp(k:(349+k));
%                         count(a+1)=count(a+1)+1;
%                     end
                             
                end
            end
        end

        for j=1:10
            subplot(5,2,j)
            hold on
            plot(filters(:,j)/count(j)*6,'r')
        end
        
        filters=zeros(500,10);
        count=zeros(1,10);
        for i=97:2:120%:trialsL
            
            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            spikes=convolved(spikes,40,62000);
            %             spikes=spikes-startTimes(i);
            spikes=spikes(121:end-120);
            bfr=mean(spikes(50:1950));
            spikes(1:startTimes(i))=[];
            
            if bfr<10
                bfr=mean(spikes);
            end
            m=bfr/5;
            m=0:m:(10*m);
            m(end)=1000;

            flickerTmp=flicker(:,i);


            for t=1:dim  
                cc=mod((t)-1,2)+1;
                if length(spikes)>10                    
                    
                    startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
                    endPoint=correctedProtocols(end,1,i)+1;
                    
                    for k=startPoint:5:endPoint-1000
                        a=mean(spikes(350+k:400+k));
                        
                        filters(:,find(m>a,1)-1)=filters(:,find(m>a,1)-1)+flickerTmp(k:(499+k));
                        count(find(m>a,1)-1)=count(find(m>a,1)-1)+1;
                    end

                             
                end
            end
        end
        for j=1:10
            subplot(5,2,j)
            hold on
            plot(filters(:,j)/count(j),'c')
        end        
        
        filters=zeros(500,10);
        count=zeros(1,10);
        for i=98:2:120%:trialsL
            
            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            spikes=convolved(spikes,40,62000);
            %             spikes=spikes-startTimes(i);
            spikes=spikes(121:end-120);
            bfr=mean(spikes(50:1950));
            spikes(1:startTimes(i))=[];
            
            if bfr<10
                bfr=mean(spikes);
            end
            m=bfr/5;
            m=0:m:(10*m);
            m(end)=1000;

            flickerTmp=flicker(:,i);


            for t=1:dim  
                cc=mod((t)-1,2)+1;
                if length(spikes)>10                    
                    
                    startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
                    endPoint=correctedProtocols(end,1,i)+1;
                    
                    for k=startPoint:5:endPoint-1000
                        a=mean(spikes(350+k:400+k));
                        
                        filters(:,find(m>a,1)-1)=filters(:,find(m>a,1)-1)+flickerTmp(k:(499+k));
                        count(find(m>a,1)-1)=count(find(m>a,1)-1)+1;
                    end

                             
                end
            end
        end
        for j=1:10
            subplot(5,2,j)
            hold on
            plot(filters(:,j)/count(j)*4,'m')
        end           
        
        
        
        
        cnt=6
        load([mainpath,'units/',units(cnt).name]);    
        filters=zeros(1000,1);
        for i=49:2:72%:trialsL
            
            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            spikes=convolved(spikes,40,62000);
            spikes=spikes(121:end-120);
            spikes(1:startTimes(i))=[];
            spikes=diff(spikes);
            flickerTmp=flicker(:,i);
            startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
            endPoint=correctedProtocols(end,1,i)+1;
            
            for k=startPoint:endPoint-1500
                a=mean(spikes(550+k));                
                filters=filters+flickerTmp(k:(999+k))*a;
            end
        end
        figure
        plot(filters)
        hold on
        filters=zeros(1000,1);
        for i=50:2:72%:trialsL
            
            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            spikes=convolved(spikes,40,62000);
            spikes=spikes(121:end-120);
            spikes(1:startTimes(i))=[];
            spikes=diff(spikes);
            flickerTmp=flicker(:,i);
            startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
            endPoint=correctedProtocols(end,1,i)+1;
            
            for k=startPoint:5:endPoint-1500
                a=mean(spikes(550+k:590+k));                
                filters=filters+flickerTmp(k:(999+k))*a;
            end
        end    
        plot(filters*6,'r')
         filters=zeros(1000,1);
        for i=97:2:120%:trialsL
            
            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            spikes=convolved(spikes,40,62000);
            spikes=spikes(121:end-120);
            spikes(1:startTimes(i))=[];
            spikes=diff(spikes);
            flickerTmp=flicker(:,i);
            startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
            endPoint=correctedProtocols(end,1,i)+1;
            
            for k=startPoint:5:endPoint-1500
                a=mean(spikes(550+k:590+k));                
                filters=filters+flickerTmp(k:(999+k))*a;
            end
        end
        plot(filters,'c')
        hold on
        filters=zeros(1000,1);
        for i=98:2:120%:trialsL
            
            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            spikes=convolved(spikes,40,62000);
            spikes=spikes(121:end-120);
            spikes(1:startTimes(i))=[];
            spikes=diff(spikes);
            flickerTmp=flicker(:,i);
            startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
            endPoint=correctedProtocols(end,1,i)+1;
            
            for k=startPoint:5:endPoint-1500
                a=mean(spikes(550+k:590+k));                
                filters=filters+flickerTmp(k:(999+k))*a;
            end
        end    
        plot(filters*6,'m')       
        
      
        
        
        
        
        
        
        
        
        
wind=435;
for cnt=1:39
    cnt=9
    load([mainpath,'units/',units(cnt).name]);
    
    bfr=0;
    for i=49:72%:trialsL
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        bfr=bfr+sum(spikes<startTimes(i))/24;
    end
    bfr1=bfr/(startTimes(i)/50);
    
    bfr=0;
    for i=97:120%:trialsL
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        bfr=bfr+sum(spikes<startTimes(i))/24;
    end
    bfr=bfr/(startTimes(i)/50);
    
    filters=zeros(500,10);
    count=zeros(1,10);
    
    for i=49:2:72%:trialsL
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        spikes=spikes-startTimes(i);
        spikes(spikes<500)=[];
        
        flickerTmp=flicker(:,i);
        
        startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
        endPoint=correctedProtocols(end,1,i)+1;
        
        for k=startPoint:17:endPoint-1000
            a=sum(spikes>=(400+k)&spikes<(wind+k));
            if a>9
                a=9;
            end
            filters(500:-1:1,a+1)=filters(500:-1:1,a+1)+flickerTmp(k:(499+k));
            count(a+1)=count(a+1)+1;
        end
    end
    figure
    for j=1:10
        subplot(5,2,j)
        hold on
        plot(filters(:,j)/count(j))
        title([int2str(cnt),'  ND6 ',int2str(bfr1),' ND4 ',int2str(bfr)])        
    end
    
    filters=zeros(500,10);
    count=zeros(1,10);
    
    for i=50:2:72%:trialsL
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        spikes=spikes-startTimes(i);
        spikes(spikes<500)=[];
        
        flickerTmp=flicker(:,i);
        
        startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
        endPoint=correctedProtocols(end,1,i)+1;
        
        
        for k=startPoint:17:endPoint-1000
            a=sum(spikes>=(400+k)&spikes<(wind+k));
            if a>9
                a=9;
            end
            filters(500:-1:1,a+1)=filters(500:-1:1,a+1)+flickerTmp(k:(499+k));
            count(a+1)=count(a+1)+1;
        end

    end
    for j=1:10
        subplot(5,2,j)
        hold on
        plot(filters(:,j)/count(j)*6,'r')
    end

    
%     filters=zeros(500,10);
%     count=zeros(1,10);
%     
%     for i=97:2:120%:trialsL
%         spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
%         spikes=spikes-startTimes(i);
%         spikes(spikes<500)=[];
%         
%         flickerTmp=flicker(:,i);
%         
%         startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
%         endPoint=correctedProtocols(end,1,i)+1;
%         
%         
%         for k=startPoint:17:endPoint-1000
%             a=sum(spikes>=(400+k)&spikes<(wind+k));
%             if a>9
%                 a=9;
%             end
%             filters(:,a+1)=filters(:,a+1)+flickerTmp(k:(499+k));
%             count(a+1)=count(a+1)+1;
%         end
% 
% 
%     end
%     for j=1:10
%         subplot(5,2,j)
%         hold on
%         plot(filters(:,j)/count(j),'c')
%     end
%     
%     
%     filters=zeros(500,10);
%     count=zeros(1,10);
%     
%     for i=98:2:120%:trialsL
%         spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
%         spikes=spikes-startTimes(i);
%         spikes(spikes<500)=[];
%         
%         flickerTmp=flicker(:,i);
%         
%         startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
%         endPoint=correctedProtocols(end,1,i)+1;
%         
%         
%         for k=startPoint:17:endPoint-1000
%             a=sum(spikes>=(400+k)&spikes<(wind+k));
%             if a>9
%                 a=9;
%             end
%             filters(:,a+1)=filters(:,a+1)+flickerTmp(k:(499+k));
%             count(a+1)=count(a+1)+1;
%         end
% 
%     end
%     for j=1:10
%         subplot(5,2,j)
%         hold on
%         plot(filters(:,j)/count(j)*6,'m')
%         line([0,500],[0,0],'color','k')
%         axis tight
%     end
    set(gcf,'position',[1400 1 440 952])
    saveas(gcf,['/mnt/muench_data/data/alexandra/MEA_data/analysis/misc_pics/',int2str(cnt),'.png'])
    close(gcf)
end





date='20121023'

load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])
load([path2save,'protocols_',codeWord])

if length(unique(startTimes))>1
    for i=1:size(correctedProtocols,3)
        correctedProtocols(:,1,i)=correctedProtocols(:,1,i)-startTimes(i);
    end
else
    correctedProtocols(:,1,:)=correctedProtocols(:,1,:)-startTimes(1);
end



figure
wind=500;
sigDif=5;
t=1;
st=73
for cnt=1:10
    load([mainpath,'units/',units(cnt).name]);
    
    bfr=0;
    hc_fr=0;
    lc_fr=0;
    for i=st:st+23%:trialsL
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        bfr=bfr+sum(spikes<startTimes(i))/24;
        if mod(i,2)>0
            hc_fr=hc_fr+sum(spikes>startTimes(i)&spikes<57536)/(12*29);
        else
            lc_fr=hc_fr+sum(spikes>startTimes(i)&spikes<57536)/(12*29);
        end
    end
       
    filters=zeros(500,10);
    count=zeros(1,10);

    for i=st:2:st+23%:trialsL
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        spikes=spikes-startTimes(i);
        spikes(spikes<500)=[];
        
        flickerTmp=flicker(:,i);
        
        startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
        endPoint=correctedProtocols(end,1,i)+1;
        
        for k=startPoint:16:endPoint-1000

%         for mm=1:3500
%             k=correctedProtocols(per*(t-1)+mm,1,i)+1;
            a=sum(spikes>=(400+k)&spikes<(wind+k));
            if a>9
                a=9;
            end
            filters(500:-1:1,a+1)=filters(500:-1:1,a+1)+flickerTmp(k:(499+k));
            count(a+1)=count(a+1)+1;
        end
    end
    mk=count;
    cc=cnt;
    while cc>10
        cc=cc-10;
    end
    for j=1:6
        subplot(7,10,cc+(j-1)*10)
        plot(filters(:,j)/count(j))
        hold on
    end
    
    filters=zeros(500,10);
    count=zeros(1,10);
    
    for i=st+1:2:st+23%:trialsL
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        spikes=spikes-startTimes(i);
        spikes(spikes<500)=[];
        
        flickerTmp=flicker(:,i);
        
        startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
        endPoint=correctedProtocols(end,1,i)+1;
        
        
        for k=startPoint:16:endPoint-1000

%         for mm=1:3500
%             k=correctedProtocols(per*(t-1)+mm,1,i)+1;
            a=sum(spikes>=(400+k)&spikes<(wind+k));
            if a>9
                a=9;
            end
            filters(500:-1:1,a+1)=filters(500:-1:1,a+1)+flickerTmp(k:(499+k));
            count(a+1)=count(a+1)+1;
        end

    end
    cc=cnt;
    while cc>10
        cc=cc-10;
    end
    for j=1:6
        subplot(7,10,cc+(j-1)*10)        
        plot(filters(:,j)/count(j)*sigDif,'r')
        line([0 500],[0 0],'color','k')
        
        title(['h ',int2str(mk(j)),'  l ',int2str(count(j))]) 
        axis tight
    end
    
    a=mean(LinearFilter(:,1,st:2:st+23,cnt),3)/mean(SpikeCount(st:2:st+23,cnt));
    b=mean(LinearFilter(:,1,st+1:2:st+23,cnt),3)/mean(SpikeCount(st+1:2:st+23,cnt));
    cc=cnt;
    while cc>10
        cc=cc-10;
    end
    subplot(7,10,cc+6*10)
    plot(a)
    hold on
    plot(b*sigDif,'r')
    line([0 500],[0 0],'color','k')
    title([int2str(cnt),' b ',int2str(bfr),' h ',int2str(hc_fr),' l ',int2str(lc_fr)])
       
    axis tight
    drawnow
end








figure
wind=500;
sigDif=5;
t=1;
st=49
col='b';
for cnt=1:39
    load([mainpath,'units/',units(cnt).name]);
    
    bfr=0;
    hc_fr=0;
    lc_fr=0;
    for i=st:st+23%:trialsL
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        bfr=bfr+sum(spikes<startTimes(i))/24;
        if mod(i,2)>0
            hc_fr=hc_fr+sum(spikes>startTimes(i)&spikes<57536)/(12*29);
        else
            lc_fr=hc_fr+sum(spikes>startTimes(i)&spikes<57536)/(12*29);
        end
    end
       
    filters=zeros(500,10);
    count=zeros(1,10);

    for i=st:2:st+23%:trialsL
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        spikes=spikes-startTimes(i);
        spikes(spikes<500)=[];
        
        flickerTmp=flicker(:,i);
        
        startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
        endPoint=correctedProtocols(end,1,i)+1;
        
        for k=startPoint:16:endPoint-1000

%         for mm=1:3500
%             k=correctedProtocols(per*(t-1)+mm,1,i)+1;
            a=sum(spikes>=(400+k)&spikes<(wind+k));
            if a>9
                a=9;
            end
            filters(500:-1:1,a+1)=filters(500:-1:1,a+1)+flickerTmp(k:(499+k));
            count(a+1)=count(a+1)+1;
        end
    end

    for j=1:10
        filters(:,j)=filters(:,j)/count(j);
    end
    
    filters_low=zeros(500,10);
    count_low=zeros(1,10);    
    for i=st+1:2:st+23%:trialsL
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        spikes=spikes-startTimes(i);
        spikes(spikes<500)=[];
        
        flickerTmp=flicker(:,i);
        
        startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
        endPoint=correctedProtocols(end,1,i)+1;
        
        
        for k=startPoint:16:endPoint-1000

%         for mm=1:3500
%             k=correctedProtocols(per*(t-1)+mm,1,i)+1;
            a=sum(spikes>=(400+k)&spikes<(wind+k));
            if a>9
                a=9;
            end
            filters_low(500:-1:1,a+1)=filters_low(500:-1:1,a+1)+flickerTmp(k:(499+k));
            count_low(a+1)=count_low(a+1)+1;
        end

    end

    for j=1:10   
        filters_low(:,j)=filters_low(:,j)/count_low(j);
    end
    
    
    
    clear a a1 b b1 c c1 f
    for j=1:10
        f(j)=corr(filters(60:300,j),filters_low(60:300,j)*5);
        if f(j)>0.6
            a=sum(abs(filters(:,j)));
            b=sum(abs(filters_low(:,j)));
            c(j)=a/b;
            for ccr=1:100
                ck=(c(j)-5)+0.1*ccr;
                tmm(ccr)=sqrt(sum((filters(60:200,j)-filters_low(60:200,j)*ck).*(filters(60:200,j)-filters_low(60:200,j)*ck)));
            end
            [c1(j) ccr]=min(tmm);
            c(j)=(c(j)-5)+0.1*ccr;
            if isnan(c(j))
                disp('NANA!')
                cnt
            end
        else
            c1(j)=nan;
            c(j)=nan;
        end
    end
    subplot(6,7,cnt)
    plot(c,col)
    hold on
    plot(c,'*','color',col)
    line([0 10],[5,5],'color','k')   
    drawnow

    
end

st=49
for cnt=1:39
    subplot(6,7,cnt)
    axis([0 11 0 40])
%     line([0 10],[1,1],'color','k') 
    load([mainpath,'units/',units(cnt).name]);
    bfr=0;
    hc_fr=0;
    lc_fr=0;
    for i=st:st+23%:trialsL
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        bfr=bfr+sum(spikes<startTimes(i))/24;
        if mod(i,2)>0
            hc_fr=hc_fr+sum(spikes>startTimes(i)&spikes<57536)/(12*29);
        else
            lc_fr=hc_fr+sum(spikes>startTimes(i)&spikes<57536)/(12*29);
        end
    end
    title([int2str(cnt),' b ',int2str(bfr/20),' h ',int2str(hc_fr/20),' l ',int2str(lc_fr/20)])
end







figure
wind=500;
sigDif=5;
t=1;
st=97
col='b';
linear_high=zeros(19,cnt);
linear_low=zeros(19,cnt);
for cnt=1:39
    load([mainpath,'units/',units(cnt).name]);
    
    bfr=0;
    hc_fr=0;
    lc_fr=0;
    spikesAll=zeros(60000,24);
    for i=st:st+23%:trialsL
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        tmp=convolved(spikes,40,63000);
        tmp=tmp(121:62120);
        spikesAll(:,i-st+1)=tmp(startTimes(i):startTimes(i)+60000-1);
        bfr=bfr+tmp(1:startTimes(i))/24;
    end
    bfr=mean(bfr);
    a=spikesAll(:,1:2:end);
    a=a(:);
    a=reshape(a,100,600*12);
    a=mean(a);
    max1=max(a);
    min1=min(a);
    a=spikesAll(:,2:2:end);
    a=a(:);
    a=reshape(a,100,600*12);
    a=mean(a);
    max2=max(max(a),max1);
    min2=min(min(a),min1);
    binSize=(max2-min2)/20;
    bins=min2:binSize:max2;
    bfrBin=bfr/binSize;
    bfrBin=ceil(bfrBin-bins(1)/binSize);
    r0(cnt)=bfrBin;
   
       
    filters=zeros(500,24,20);
    count=zeros(24,20);
    
    startPoint=correctedProtocols(per*(t-1)+1,1,st)+1;
    endPoint=correctedProtocols(end,1,st)+1;
    
    for k=startPoint:16:endPoint-1000
        a=mean(spikesAll(400+k:500+k,:));
        a=a/binSize;
        a=ceil(a-bins(1)/binSize);
        a(a>20)=20;
        a(a<1)=1;
        for p=1:24
            filters(500:-1:1,p,a(p))=filters(500:-1:1,p,a(p))+flicker(k:(499+k),st+p-1);
            count(p,a(p))=count(p,a(p))+1;
        end
    end

    for j=1:20
        for k=1:24
            filters(:,k,j)=filters(:,k,j)/count(k,j);
        end
    end
    filt_high=zeros(500,20);
    filt_low=filt_high;
    for str=1:20
        filt_high(:,str)=(nanmean(filters(:,1:2:end,str),2));
        filt_low(:,str)=(nanmean(filters(:,2:2:end,str),2));
    end
%     
%     
%     figure
%     for str=1:20
%         subplot(4,5,str)
%         plot(nanmean(filters(:,1:2:end,str),2))
%         hold on
%         plot(nanmean(filters(:,2:2:end,str),2)*5,'r')
%     end
%     

    % check question 1
    
    for j=1:19
        if j<bfrBin
            linear_high(j,cnt)=sum(abs(filt_high(1:300,j)))-sum(abs(filt_high(1:300,j+1)));
            linear_low(j,cnt)=sum(abs(filt_low(1:300,j)))-sum(abs(filt_low(1:300,j+1)));
        else
            linear_high(j,cnt)=sum(abs(filt_high(1:300,j+1)))-sum(abs(filt_high(1:300,j)));
            linear_low(j,cnt)=sum(abs(filt_low(1:300,j+1)))-sum(abs(filt_low(1:300,j)));
        end
%         linear_high(j,cnt)=sqrt(sum((filt_high(1:300,j+1)-filt_high(1:300,j)).*(filt_high(1:300,j+1)-filt_high(1:300,j))));
%         linear_low(j,cnt)=sqrt(sum((filt_low(1:300,j+1)-filt_low(1:300,j)).*(filt_low(1:300,j+1)-filt_low(1:300,j))));
    end

    % check question 2
    clear a a1 b b1 c c1 f
    for j=1:20
        f(j)=corr(filt_high(1:300,j),filt_low(1:300,j));
        if f(j)>0.6
            a=sum(abs(filt_high(:,j)));
            b=sum(abs(filt_low(:,j)));
            c(j)=a/b;
            for ccr=1:100
                ck=(c(j)-5)+0.1*ccr;
                tmm(ccr)=sqrt(sum((filt_high(1:300,j)-filt_low(1:300,j)*ck).*(filt_high(1:300,j)-filt_low(1:300,j)*ck)));
            end
            [c1(j) ccr]=min(tmm);
            c(j)=(c(j)-5)+0.1*ccr;
            if isnan(c(j))
                disp('NAN!')
                cnt
            end
        else
            c1(j)=nan;
            c(j)=nan;
        end
    end
%     subplot(6,7,cnt)
%     plot(c,col)
%     hold on
%     plot(c,'*','color',col)
%     line([bfrBin,bfrBin],[0,25],'color',col)
%     line([0 21],[5,5],'color','k')
%     axis([0 21 0 25])
%     drawnow    
%     title([int2str(cnt),' b ',int2str(bfr),' min ',int2str(min2),' max ',int2str(max2)])
end

%fig U1
figure
for cnt=1:39
    subplot(6,7,cnt)
    plot(linear_high(:,cnt),'-*')
    hold on
    
    plot(linear_low(:,cnt)*5,'-*r')
    axis tight
    line([r0(cnt),r0(cnt)],[-20,30],'color','k')
    
    title(int2str(cnt))
end

