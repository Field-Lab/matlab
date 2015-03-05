cd('/mnt/muench_data/user/alexandra/scripts')

% prepare matrix of firing rate and flickers

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
    metastat=zeros(trialsL+trialsL*ifDim,2); % ND, 0 - high, 1 - low contrast
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
    
    
    FR=zeros(maxLength,min(size(metastat,1),192),length(units));
    ST=FR;
    basalFR=zeros(min(size(metastat,1),192),length(units));
    basalFRstd=basalFR;
    maxLast=1;
    for cnt=1:length(units)
        cnt
        load([mainpath,'units/',units(cnt).name]);
        dimCnt=1;
        for i=1:trialsL

            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            convSpikes=convolved(spikes,40,maxLength);
            convSpikes=convSpikes(121:end-120);
            spontMean=mean(convSpikes(50:startTimes(i)-50));
            spontStd=std(convSpikes(50:startTimes(i)-50));
            convSpikes=convSpikes(startTimes(i):end);
            flickerTmp=flicker(:,i);
            flickers=zeros(maxLength,ifDim+1);
            firingRates=zeros(maxLength,ifDim+1);
            last=ones(1,ifDim+1);

            for t=1:dim  
                cc=mod((t)-1,2)+1;
                if length(spikes)>10
                    startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
                    endPoint=min(correctedProtocols(per*(t-1)+dur,1,i),length(convSpikes)); 
                    tmp=flickerTmp(startPoint:endPoint);
                    flickers(last(cc):last(cc)+length(tmp)-1,cc)=tmp;
                    tmp=convSpikes(startPoint:endPoint);
                    firingRates(last(cc):last(cc)+length(tmp)-1,cc)=tmp;
                    last(cc)=last(cc)+length(tmp);
                end
            end
            firingRates=firingRates(1:max(last)-1,:);
            flickers=flickers(1:max(last)-1,:);
            
            FR(1:max(last)-1,dimCnt,cnt)=firingRates(1:max(last)-1,1);
            ST(1:max(last)-1,dimCnt,cnt)=flickers(1:max(last)-1,1);
            
            basalFR(dimCnt,cnt)=spontMean;
            basalFRstd(dimCnt,cnt)=spontStd;
            
            if ifDim>0
                FR(1:max(last)-1,dimCnt+1,cnt)=firingRates(1:max(last)-1,2);
                ST(1:max(last)-1,dimCnt+1,cnt)=flickers(1:max(last)-1,2);                
                basalFR(dimCnt+1,cnt)=spontMean;
                basalFRstd(dimCnt+1,cnt)=spontStd;
            end
            maxLast=max(maxLast,max(last)-1);
            dimCnt=dimCnt+1+ifDim*1;
        end        
    end
    FR=FR(1:maxLast,:,:);
    ST=ST(1:maxLast,:,:);
    
    save([path2save,'fainGrainMatrices.mat'],'FR','ST','basalFR','basalFRstd','metastat','names','onOff')
    toc
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

t=1;
st=97
col='b';
linear_high=zeros(19,39);
linear_low=zeros(19,39);
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
  

    % check question 1
    
    for j=1:19
        if j<bfrBin
            linear_high(j,cnt)=sum(abs(filt_high(1:300,j)))-sum(abs(filt_high(1:300,j+1)));
            linear_low(j,cnt)=sum(abs(filt_low(1:300,j)))-sum(abs(filt_low(1:300,j+1)));
        else
            linear_high(j,cnt)=sum(abs(filt_high(1:300,j+1)))-sum(abs(filt_high(1:300,j)));
            linear_low(j,cnt)=sum(abs(filt_low(1:300,j+1)))-sum(abs(filt_low(1:300,j)));
        end
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



clear
 


date='20121023'
allFilters=zeros(500,192,30,39);
allCounts=zeros(192,30,39);
for cnt=1:size(FR,3)
    cnt
    % make trials list per nd
    for i=1:8
        tmp=find(metastat(:,1)==(9-i));
        a=find(diff(tmp)>1);   
        if ~isempty(a)
            trialsPerND(:,i)=tmp(1:a);
        else
            trialsPerND(:,i)=tmp;
        end
    end
    
    % calculate bins for every trial
    
    bins=zeros(30,numel(trialsPerND));
    for j=1:numel(trialsPerND)
        binSize(j)=max(basalFRstd(j),4.9275);
        a=basalFR(j):binSize(j):(basalFR(j)+binSize(j)*15);
        b=basalFR(j)-binSize(j):-binSize(j):-(basalFR(j)+binSize(j)*15);
        binsTmp=sort([a b]);
        binsTmp(binsTmp<0)=[];
        bins(1:length(binsTmp),j)=binsTmp;
        maxBins(j)=length(binsTmp);
        bfrBin=basalFR(j)/binSize(j);
        bfrBin=ceil(bfrBin-binsTmp(1)/binSize(j));
        r0(j,cnt)=bfrBin;
    end
    bins=bins(1:max(maxBins),:);
    

    
    filters=zeros(500,numel(trialsPerND),max(maxBins));
    count=zeros(numel(trialsPerND),max(maxBins));
    
    stopPoint=find(ST(:,i,cnt)==0,1,'last');
    if ~isempty(stopPoint)
        stopPoint=size(ST,1);
    end
    
    for k=1:10:stopPoint-900;
        a=mean(FR(400+k:480+k,:,cnt));
        a=a./binSize;
        a=ceil(a-bins(1,:)./binSize);
        if sum(a>maxBins)>0
            a(a>maxBins)=maxBins(a>maxBins);
        end
        a(a<1)=1;
        for p=1:numel(trialsPerND)
            filters(500:-1:1,p,a(p))=filters(500:-1:1,p,a(p))+ST(k:(499+k),p,cnt);
            count(p,a(p))=count(p,a(p))+1;
        end
    end
    
    for j=1:max(maxBins)
        for k=1:numel(trialsPerND)
            filters(:,k,j)=filters(:,k,j)/count(k,j);
        end
    end
     
    allFilters(:,:,1:max(maxBins),cnt)=filters;
    allCounts(:,1:max(maxBins),cnt)=count;
    
end

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/fineGrainResults_',date],'allFilters','allCounts','names','onOff','r0','metastat')



clear a a1 b b1 c c1 f
    cc=1;
    for i=1:2:192
        for j=1:28
            f(j)=corr(filters(1:300,i,j),filters(1:300,i+1,j));
            if f(j)>0.6
                a=sum(abs(filters(1:300,i,j)));
                b=sum(abs(filters(1:300,i+1,j)));
                c(j)=a/b;
                for ccr=1:100
                    ck=(c(j)-5)+0.1*ccr;
                    tmm(ccr)=sqrt(sum((filters(1:300,i,j)-filters(1:300,i+1,j)*ck).*(filters(1:300,i,j)-filters(1:300,i+1,j)*ck)));
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
        coeffs(cc,:)=c;
        cc=cc+1;
    end
plot(coeffs(25:36,1:20)')
