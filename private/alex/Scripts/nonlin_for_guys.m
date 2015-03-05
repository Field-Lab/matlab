%% Nonlinearity

%% 20120329
clear
boundaries=100;
date='2012-01-30';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/!WT_Hartwig/',date,'/LF_FFFlicker.mat'])
load(['/mnt/muench_data/data/Hartwig/MEA_data/new_params_2013_04_18/linear_filters_manual.mat'])

for i=1:18
    tmp=reshape(mean(LinearFilter(:,1:2:10,i,:),2),500,size(names,2));
    subplot(6,3,i)
    plot(tmp)
end
ourExp=LF_manual(:,149:149+size(LinearFilter,4)-1);
% nds=[8 8 8 8 8 8 6 6 6 6 4 4 4 4 3 3 3 3];
nds=[1 1 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];


%% Prepare nonlin (all trials and dimensions)
cd('/mnt/muench_data/user/alexandra/scripts')
clear

dates=cell(15,1);
dates{1}='2012-01-30';

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

rangeHC=-0.3:0.3/30:0.3;
rangeLC=-0.07:0.07/30:0.07;

for datesCNT=1:15
    date=dates{datesCNT}
    
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/!WT_Hartwig/',date,'/LF_FFFlicker.mat'])
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/!WT_Hartwig/',date,'/protocols_FFFlicker.mat'])
        
    load([path2save,date,'_nonlinear'],'nonLinear','names')
    nonlin=zeros(61,dim,size(LinearFilter,3),length(units));    
    
    convSpikes=zeros(size(flicker,1),size(LinearFilter,3),length(units));
    for cnt=1:length(units)
        load([mainpath,'units/',units(cnt).name]);
        for i=1:size(LinearFilter,3)
            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            spikes=spikes-startTimes(i)+1;
            spikes(spikes<500|spikes>size(flicker,1))=[];
            tmp=convolved(spikes,40,size(flicker,1));
            convSpikes(:,i,cnt)=tmp(121:end-120);
        end
    end    
    
    for i=1:size(LinearFilter,3)
        i
        for t=1:dim
            subflicker=flicker(correctedProtocols(per*(t-1)+1,1,i)-startTimes(i)+1:correctedProtocols(per*(t-1)+dur,1,i)-startTimes(i),i);
            tmp=zeros(length(subflicker)-499,length(units));
            a=reshape(LinearFilter(500:-1:1,t,i,:),500,length(units))';
            for kk=1:size(a,1)
                a(kk,:)=a(kk,:)-mean(a(kk,:));
                a(kk,:)=a(kk,:)./sum(abs(a(kk,:)));
            end
            for rt=500:length(subflicker)
                tmp(rt-499,:)=a*subflicker(rt-499:rt);
            end
            
            a=convSpikes((correctedProtocols(per*(t-1)+1,1,i)-startTimes(i))+500:(correctedProtocols(per*(t-1)+dur,1,i)-startTimes(i)),i,:);
            a=reshape(a,size(a,1),size(a,3));
            if std(subflicker)>0.15
                
                x=rangeHC;
                stepSize=0.3/30;
            else
                x=rangeLC;
                stepSize=0.07/30;
            end
            
            kk=1;
            for rt=x
                tt=tmp>=rt&tmp<rt+stepSize;
                [timeP,cellID]=find(tt);
                for l=unique(cellID)'
                    nonlin(kk,t,i,l)=mean(a(timeP(cellID==l),l));
                end
                kk=kk+1;
            end
        end
    end
    save([path2save,date,'_nonlin'],'nonlin','names','rangeHC','rangeLC')
end


%% reduce dimensions and trials, smooth out noise

cd('/mnt/muench_data/user/alexandra/scripts')
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
    date=dates{datesCNT}
    path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
    load([path2save,date,'_nonlin'],'nonlin','names','rangeHC','rangeLC')
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date],'units')

    switch size(nonlin,3)
        case 192
            
            nonl=zeros(61,64,length(units));
            cc=1;
            nonlin=reshape(nonlin,61,192,size(nonlin,4));
            nonlin(nonlin==0)=nan;
            for i=1:6:192
                nonl(:,cc,:)=nanmean(nonlin(:,i:2:i+5,:),2);
                nonl(:,cc+1,:)=nanmean(nonlin(:,i+1:2:i+5,:),2);
                cc=cc+2;
            end
            form=zeros(1,192); % low contrast
            form(1:2:end)=1; % high contrast
        case 46
            nonl=zeros(61,32,length(units));
            cc=1;
            nonlin(nonlin==0)=nan;
            tmpHC=nanmean(nonlin(:,1:2:end,:,:),2);
            tmpLC=nanmean(nonlin(:,2:2:end,:,:),2);
            tmpHC=reshape(tmpHC,61,size(nonlin,3),size(nonlin,4));
            tmpLC=reshape(tmpLC,61,size(nonlin,3),size(nonlin,4));
            for i=1:2:32
                nonl(:,cc,:)=nanmean(tmpHC(:,i:i+1,:),2);
                nonl(:,cc+1,:)=nanmean(tmpLC(:,i:i+1,:),2);
                cc=cc+2;
            end
            form=zeros(1,size(nonl,2)); % low contrast
            form(1:2:end)=1; % high contrast
        case 82
            strt=10;
            if strcmp(date,'20130220')
                strt=11;
            end            
            nonl=zeros(61,24,length(units));
            cc=1;
            nonlin=reshape(nonlin,61,size(nonlin,3),size(nonlin,4));
            nonlin(nonlin==0)=nan;
            for i=strt:3:size(nonlin,2)-1
                nonl(:,cc,:)=nanmean(nonlin(:,i:i+2,:),2);
                cc=cc+1;
            end
            form=ones(1,size(nonl,2)); % high contrast
        case 81
            nonl=zeros(61,24,length(units));
            cc=1;
            nonlin=reshape(nonlin,61,size(nonlin,3),size(nonlin,4));
            nonlin(nonlin==0)=nan;
            for i=10:3:size(nonlin,2)
                nonl(:,cc,:)=nanmean(nonlin(:,i:i+2,:),2);                
                cc=cc+1;
            end
            form=ones(1,size(nonl,2)); % high contrast

        case 16
            nonl=zeros(61,32,length(units));
            cc=1;
            nonlin(nonlin==0)=nan;
            tmpHC=nanmean(nonlin(:,1:2:end,:,:),2);
            tmpLC=nanmean(nonlin(:,2:2:end,:,:),2);
            tmpHC=reshape(tmpHC,61,size(nonlin,3),size(nonlin,4));
            tmpLC=reshape(tmpLC,61,size(nonlin,3),size(nonlin,4));
            for i=1:16
                nonl(:,cc,:)=tmpHC(:,i,:);
                nonl(:,cc+1,:)=tmpLC(:,i,:);
                cc=cc+2;
            end
            form=zeros(1,size(nonl,2)); % low contrast
            form(1:2:end)=1; % high contrast
        case 24
            nonl=zeros(61,48,length(units));
            cc=1;
            nonlin(nonlin==0)=nan;
            tmpHC=nanmean(nonlin(:,1:2:end,:,:),2);
            tmpLC=nanmean(nonlin(:,2:2:end,:,:),2);
            tmpHC=reshape(tmpHC,61,size(nonlin,3),size(nonlin,4));
            tmpLC=reshape(tmpLC,61,size(nonlin,3),size(nonlin,4));
            for i=1:24
                nonl(:,cc,:)=tmpHC(:,i,:);
                nonl(:,cc+1,:)=tmpLC(:,i,:);
                cc=cc+2;
            end
            form=zeros(1,size(nonl,2)); % low contrast
            form(1:2:end)=1; % high contrast
        case 49
            nonl=zeros(61,32,length(units));
            cc=1;
            nonlin(nonlin==0)=nan;
            tmpHC=nanmean(nonlin(:,1:2:end,:,:),2);
            tmpLC=nanmean(nonlin(:,2:2:end,:,:),2);
            tmpHC=reshape(tmpHC,61,size(nonlin,3),size(nonlin,4));
            tmpLC=reshape(tmpLC,61,size(nonlin,3),size(nonlin,4));
            for i=1:3:48
                nonl(:,cc,:)=nanmean(tmpHC(:,i:i+2,:),2);
                nonl(:,cc+1,:)=nanmean(tmpLC(:,i:i+2,:),2);
                cc=cc+2;
            end
            form=zeros(1,size(nonl,2)); % low contrast
            form(1:2:end)=1; % high contrast
        case 42
            nonl=zeros(61,28,length(units));
            cc=1;
            nonlin(nonlin==0)=nan;
            tmpHC=nanmean(nonlin(:,1:2:end,:,:),2);
            tmpLC=nanmean(nonlin(:,2:2:end,:,:),2);
            tmpHC=reshape(tmpHC,61,size(nonlin,3),size(nonlin,4));
            tmpLC=reshape(tmpLC,61,size(nonlin,3),size(nonlin,4));
            for i=1:3:42
                nonl(:,cc,:)=nanmean(tmpHC(:,i:i+2,:),2);
                nonl(:,cc+1,:)=nanmean(tmpLC(:,i:i+2,:),2);
                cc=cc+2;
            end
            form=zeros(1,size(nonl,2)); % low contrast
            form(1:2:end)=1; % high contrast
        case 36
            nonl=zeros(61,32,length(units));
            cc=1;
            nonlin(nonlin==0)=nan;
            tmpHC=nanmean(nonlin(:,1:2:end,:,:),2);
            tmpLC=nanmean(nonlin(:,2:2:end,:,:),2);
            tmpHC=reshape(tmpHC,61,size(nonlin,3),size(nonlin,4));
            tmpLC=reshape(tmpLC,61,size(nonlin,3),size(nonlin,4));
            for i=5:2:36
                nonl(:,cc,:)=nanmean(tmpHC(:,i:i+1,:),2);
                nonl(:,cc+1,:)=nanmean(tmpLC(:,i:i+1,:),2);
                cc=cc+2;
            end
            form=zeros(1,size(nonl,2)); % low contrast
            form(1:2:end)=1; % high contrast
    end
    save([path2save,date,'_nonl'],'nonl','names','rangeHC','rangeLC', 'form')    
end




%% do fitting

cd('/mnt/muench_data/user/alexandra/scripts')
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

g = fittype('a+b/(1+exp((c-x)/d))');

for datesCNT=1:5
    date=dates{datesCNT}
    path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];    
    load([path2save,date,'_nonl'])
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date],'names')
    nonLinear=zeros(4,size(nonl,2),size(nonl,3));
    
    for i=1:size(nonl,2)
        i
        if form(i)
            x=rangeHC;
        else
            x=rangeLC;
        end
        for cnt=1:size(nonl,3)
            a=nonl(:,i,cnt);
            if sum(isnan(a))<30
                b=find(isnan(a(1:30)));
                if ~isempty(b)
                    a(b)=min(a(b(end)+1:b(end)+4));
                end
                b=find(isnan(a(30:end)))+29;
                if ~isempty(b)
                    c=sort(a(~isnan(a)));
                    c=mean(c(end-5:end));
                    a(b)=c;
                end
                
                try
                    tmp_f=max(a)-min(a);
                    tmp_f=find(a>tmp_f,1);
                    res=fit(x',a,g,'StartPoint', [min(a),max(a),x(tmp_f),0.03]);
                    nonLinear(1,i,cnt)=res.a;
                    nonLinear(2,i,cnt)=res.b;
                    nonLinear(3,i,cnt)=res.c;
                    nonLinear(4,i,cnt)=res.d;
                catch
                    disp('could NOT fit')
                    date
                    cnt
                    disp('END OF could NOT fit')
                end
            end
        end
    end
    save([path2save,'gain'],'nonLinear','form','names')
end


