
% strategies

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
    
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date],'per','mainpath','path2save','codeWord','dim','dur','units')
    load([path2save,'protocols_',codeWord])   
    load([path2save,'convFilters_std'],'posFilter_std','negFilter_std','linFilter_noDim','basalFR','basalFRstd','metastat')
    
    
    if length(unique(startTimes))>1
        for i=1:size(correctedProtocols,3)
            correctedProtocols(:,1,i)=correctedProtocols(:,1,i)-startTimes(i);
        end
    else
        correctedProtocols(:,1,:)=correctedProtocols(:,1,:)-startTimes(1);
    end
    
    trialsL=min(length(file_list),192);
    
    
    ifDim=double(dim>1);
    binSize=50;
    a=floor(maxLength/binSize);
    m=zeros(1,a*binSize);
    binSp=zeros(a,trialsL,cnt);
    for cnt=1:length(units)
        cnt
        load([mainpath,'units/',units(cnt).name]);
        for i=1:trialsL
            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            m(spikes(spikes<a*binSize))=1;
            binSp(:,i,cnt)=sum(reshape(m,binSize,a));    
            m(m==1)=0;
        end
    end
    
    if ifDim>0
        
    end
    
    bar(binSp(:,46,5))
             
%     save([path2save,'convFilters_1std'],'posFilter_std','negFilter_std','linFilter_noDim','basalFR','basalFRstd','metastat')
    toc
end