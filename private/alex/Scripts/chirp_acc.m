%% Accumulate
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
codeWord='chirp';
chirp=zeros(23000,32,520);
names=cell(520,1);
real_chirp=zeros(520,1);
cellID=1;
for datesCNT=1:15
    date=dates{datesCNT}
    if strcmp(date,'20120714')
        bord=1:35;
        st=4;
        bord(1:5:35)=[];
    elseif datesCNT>5&datesCNT<10
        bord=sort([5:8:64 6:8:64 7:8:64 8:8:64]);
        st=0;
    elseif datesCNT>9&datesCNT<12
        bord=sort([9:12:96 10:12:96 11:12:96 12:12:96]);
        st=0;
    else
        bord=1:40;
        bord(1:5:40)=[];
        st=0;
    end
    mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
    hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
    heka=dir([hekapath,'*.phys']);
    
   
    % select heka files
    file_list=[];
    ff=[];
    for i=1:length(heka)
        if ~isempty(regexp(heka(i).name,codeWord, 'once'))
            file_list=[file_list i];
        end
        if ~isempty(regexp(heka(i).name,'nd_8', 'once'))
            ff=[ff i];
        end
    end
    units=dir([mainpath,'units/*.mat']);
    if ~isempty(file_list)
        if ~isempty(ff)
            file_list(file_list<ff(1))=[];
        end
        
        % get protocol
        protocol=read_header_field_heka(hekapath, heka(file_list(1)).name, 'Stimulus Protocol');
        protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
        protocol=round(protocol(2:end-2,[1,4]));

        for cnt=1:length(units)
            
            load([mainpath,'units/',units(cnt).name]);
            cr=1;
            for i=bord
                spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
                conv=convolved(spikes,40,23000);
                chirp(:,cr+st,cellID)=conv(121:end-120);
                cr=cr+1;
            end
            
            if isempty(regexp(units(cnt).name,date, 'once'))
                names{cellID}=[units(cnt).name(1),date,units(cnt).name(10:end-4)];
            else
                names{cellID}=units(cnt).name(1:end-4);
            end
            real_chirp(cellID)=1;
            cellID=cellID+1;            
        end
    else
        cellID=cellID+length(units);
    end

end


% get protocol
protocol=read_header_field_heka(hekapath, heka(file_list(1)).name, 'Stimulus Protocol');
protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
protocol=protocol(2:end-2,[1,4]);
tmp=(0:size(protocol,1)-1)'*(7.85/3600);
protocol(:,1)=round(protocol(:,1)+tmp);

flips=protocol(:,1);
pxls=protocol(2:end,2);
    
flips=flips-flips(1)+1;

tmp=zeros(flips(end),1);
tmp(flips(2:end),1)=pxls;
    
droppedMax=max(diff(flips));
for j=1:droppedMax+1
    subFrame=flips(2:end)-j;
    subFrame(subFrame<1)=[];
    subFrame(find(tmp(subFrame,1)))=[];
    tmp(subFrame,1)=tmp(subFrame+1,1);
end
flicker=[ones(protocol(1,1),1)+29; tmp];
figure
plot(flicker)
hold on
[a,b]=max(flicker(11660:12110))
b=b+11660+52
plot(b:501:19523,flicker(b:501:19523),'*r')
plot(b-250:501:19523,flicker(b-250:501:19523),'*r')

maxpos=b:501:19523;
minneg=b-250:501:19523;
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')

%% Plots

for i=1:8
    subplot(4,2,i)
%     plot(protocol(:,1),protocol(:,2),'k')
    hold on
    plot(mean(mean(chirp(:,(i-1)*4+1:i*4,onOff>0&real_chirp'>0),2),3),'r')
    plot(mean(mean(chirp(:,(i-1)*4+1:i*4,onOff<0&real_chirp'>0),2),3))
    
end


figure
nds='87654321'
k=onOff>0;
for i=1:8
    subplot(4,2,i)
    hold on
    wrongCells=(frMean_HC(:,i)-frMean_spont(:,i))<-3;
    ambiCells=(frMean_HC(:,i)-frMean_spont(:,i))>=-3&(frMean_HC(:,i)-frMean_spont(:,i))<=3;
    regCells=(frMean_HC(:,i)-frMean_spont(:,i))>3;
    a=mean(chirp(:,(i-1)*4+1:i*4,k&real_chirp'>0&wrongCells'>0),2);
    a=reshape(a,23000,size(a,3));
    a=a-repmat(mean(a(500:2500,:)),23000,1);    
    a=mean(a,2);
    [ind,lags]=xcorr(a(12000:19500),flicker(12000:19500),'coeff');
    corCoef(i,1)=corr(flicker(3000:6000),a(3000:6000));
    corCoef(i,2)=corr(flicker(8000:10000),a(8000:10000));
    plot(a,'r')
%     plot(a(1:19500)-(flicker(1:19500)-30),'r')    
    for indd=1:length(maxpos)
        [maxx(indd,i,1) indmax(indd,i,1)]=max(a(maxpos(indd)-100:maxpos(indd)+100));
        indmax(indd,i,1)=indmax(indd,i,1)+maxpos(indd)-100-1;
        [minn(indd,i,1) indmin(indd,i,1)]=min(a(minneg(indd)-100:minneg(indd)+100));
        indmin(indd,i,1)=indmin(indd,i,1)+minneg(indd)-100-1;
    end
%     plot(indmax(:,i,1),maxx(:,i,1),'*m')
%     plot(indmin(:,i,1),minn(:,i,1),'*m')
    line([0 23000],[mean(a(200:2700)),mean(a(200:2700))],'color','k')
    
    a=mean(chirp(:,(i-1)*4+1:i*4,k&real_chirp'>0&regCells'>0),2);
    a=reshape(a,23000,size(a,3));
    a=a-repmat(mean(a(500:2500,:)),23000,1);   
    a=mean(a,2);
        [ind,lags]=xcorr(a(12000:19500),flicker(12000:19500),'coeff');
    corCoef(i,3)=corr(flicker(3000:6000),a(3000:6000));
    corCoef(i,4)=corr(flicker(8000:10000),a(8000:10000));
    plot(a)
    for indd=1:length(maxpos)
        [maxx(indd,i,2) indmax(indd,i,2)]=max(a(maxpos(indd)-100:maxpos(indd)+100));
        indmax(indd,i,2)=indmax(indd,i,2)+maxpos(indd)-100-1;
        [minn(indd,i,2) indmin(indd,i,2)]=min(a(minneg(indd)-100:minneg(indd)+100));
        indmin(indd,i,2)=indmin(indd,i,2)+minneg(indd)-100-1;
    end
%     plot(indmax(:,i,2),maxx(:,i,2),'*c')
%     plot(indmin(:,i,2),minn(:,i,2),'*c')
%     plot(a(1:19500)-(flicker(1:19500)-30))
    line([0 23000],[mean(a(200:2700)),mean(a(200:2700))],'color','k')
    
%     plot((flicker(1:19500)-30),'color','k');
    
%     a=mean(chirp(:,(i-1)*4+1:i*4,k&real_chirp'>0&ambiCells'>0),2);
%     a=reshape(a,23000,size(a,3));
% %     a=a-repmat(mean(a(500:2500,:)),23000,1);   
%     a=mean(a,2);
%     plot(a,'g')
%     line([0 23000],[mean(a(200:2700)),mean(a(200:2700))],'color','k')
%     
    a=mean(chirp(:,(i-1)*4+1:i*4,onOff<0&real_chirp'>0),2);
    a=reshape(a,23000,size(a,3));
    a=a-repmat(mean(a(500:2500,:)),23000,1);   
    a=mean(a,2);
    plot(a,'g') 
    for indd=1:length(maxpos)
        [minn(indd,i,3) indmin(indd,i,3)]=min(a(maxpos(indd)-200:maxpos(indd)+200));
        indmin(indd,i,3)=indmin(indd,i,3)+maxpos(indd)-200-1;
        [maxx(indd,i,3) indmax(indd,i,3)]=max(a(minneg(indd)-200:minneg(indd)+200));
        indmax(indd,i,3)=indmax(indd,i,3)+minneg(indd)-200-1;
    end
%     plot(indmax(:,i,3),maxx(:,i,3),'*k')
%     plot(indmin(:,i,3),minn(:,i,3),'*k')
    title(['ND',nds(i),'   wrong cells ',int2str(sum(k&real_chirp'>0&wrongCells'>0)),'   regular cells ',int2str(sum(k&real_chirp'>0&regCells'>0)),'   ambi cells ',int2str(sum(k&real_chirp'>0&ambiCells'>0))])
end



%contrast response function
figure
for i=1:8
    subplot(4,2,i)
    plot([flicker(minneg(end:-1:1))-30; flicker(maxpos)-30],[minn(end:-1:1,i,1); maxx(:,i,1)],'r')
    hold on
    plot([flicker(minneg(end:-1:1))-30; flicker(maxpos)-30],[minn(end:-1:1,i,2); maxx(:,i,2)])
    plot([flicker(minneg(end:-1:1))-30; flicker(maxpos)-30],[maxx(end:-1:1,i,3); minn(:,i,3)],'g')
    legend({'wrong cells','reg cells','off'},'location','southeast')
    line([-30 30],[0,0],'color','k')
    line([0,0],[-30 30],'color','k')
    title(nds(i))
%     axis([-10 10 -20 20])
end





figure
nds='87654321'
k=onOff<0;
for i=1:8
    subplot(4,2,i)
    hold on
    wrongCells=frMean_spont(:,i)>20;
    ambiCells=(frMean_spont(:,i))>=2&(frMean_spont(:,i))<=20;
    regCells=frMean_spont(:,i)<2;
    a=mean(chirp(:,(i-1)*4+1:i*4,k&real_chirp'>0&wrongCells'>0),2);
    a=reshape(a,23000,size(a,3));
    a=a-repmat(mean(a(500:2500,:)),23000,1);    
    a=mean(a,2);
    [ind,lags]=xcorr(a(12000:19500),flicker(12000:19500),'coeff');
    corCoef(i,1)=corr(flicker(3000:6000),a(3000:6000));
    corCoef(i,2)=corr(flicker(8000:10000),a(8000:10000));
    plot(a,'r')
%     plot(a(1:19500)-(flicker(1:19500)-30),'r')    
    for indd=1:length(maxpos)
        [maxx(indd,i,1) indmax(indd,i,1)]=max(a(maxpos(indd)-100:maxpos(indd)+100));
        indmax(indd,i,1)=indmax(indd,i,1)+maxpos(indd)-100-1;
        [minn(indd,i,1) indmin(indd,i,1)]=min(a(minneg(indd)-100:minneg(indd)+100));
        indmin(indd,i,1)=indmin(indd,i,1)+minneg(indd)-100-1;
    end
%     plot(indmax(:,i,1),maxx(:,i,1),'*m')
%     plot(indmin(:,i,1),minn(:,i,1),'*m')
    line([0 23000],[mean(a(200:2700)),mean(a(200:2700))],'color','k')
    
    a=mean(chirp(:,(i-1)*4+1:i*4,k&real_chirp'>0&regCells'>0),2);
    a=reshape(a,23000,size(a,3));
    a=a-repmat(mean(a(500:2500,:)),23000,1);   
    a=mean(a,2);
        [ind,lags]=xcorr(a(12000:19500),flicker(12000:19500),'coeff');
    corCoef(i,3)=corr(flicker(3000:6000),a(3000:6000));
    corCoef(i,4)=corr(flicker(8000:10000),a(8000:10000));
    plot(a)
    for indd=1:length(maxpos)
        [maxx(indd,i,2) indmax(indd,i,2)]=max(a(maxpos(indd)-100:maxpos(indd)+100));
        indmax(indd,i,2)=indmax(indd,i,2)+maxpos(indd)-100-1;
        [minn(indd,i,2) indmin(indd,i,2)]=min(a(minneg(indd)-100:minneg(indd)+100));
        indmin(indd,i,2)=indmin(indd,i,2)+minneg(indd)-100-1;
    end
%     plot(indmax(:,i,2),maxx(:,i,2),'*c')
%     plot(indmin(:,i,2),minn(:,i,2),'*c')
%     plot(a(1:19500)-(flicker(1:19500)-30))
    line([0 23000],[mean(a(200:2700)),mean(a(200:2700))],'color','k')
    
%     plot((flicker(1:19500)-30),'color','k');
    
%     a=mean(chirp(:,(i-1)*4+1:i*4,k&real_chirp'>0&ambiCells'>0),2);
%     a=reshape(a,23000,size(a,3));
% %     a=a-repmat(mean(a(500:2500,:)),23000,1);   
%     a=mean(a,2);
%     plot(a,'g')
%     line([0 23000],[mean(a(200:2700)),mean(a(200:2700))],'color','k')
%     
    a=mean(chirp(:,(i-1)*4+1:i*4,onOff<0&real_chirp'>0),2);
    a=reshape(a,23000,size(a,3));
    a=a-repmat(mean(a(500:2500,:)),23000,1);   
    a=mean(a,2);
    plot(a,'g') 
    for indd=1:length(maxpos)
        [minn(indd,i,3) indmin(indd,i,3)]=min(a(maxpos(indd)-200:maxpos(indd)+200));
        indmin(indd,i,3)=indmin(indd,i,3)+maxpos(indd)-200-1;
        [maxx(indd,i,3) indmax(indd,i,3)]=max(a(minneg(indd)-200:minneg(indd)+200));
        indmax(indd,i,3)=indmax(indd,i,3)+minneg(indd)-200-1;
    end
%     plot(indmax(:,i,3),maxx(:,i,3),'*k')
%     plot(indmin(:,i,3),minn(:,i,3),'*k')
    title(['ND',nds(i),'   wrong cells ',int2str(sum(k&real_chirp'>0&wrongCells'>0)),'   regular cells ',int2str(sum(k&real_chirp'>0&regCells'>0)),'   ambi cells ',int2str(sum(k&real_chirp'>0&ambiCells'>0))])
end




figure
nds='87654321'
col='rgb';
for i=1:8
    subplot(2,4,i)
    wrongCells=(frMean_HC(onOff>0,i)-frMean_spont(onOff>0,i))<-3;
    b=frMean_spont(onOff>0,i);
    b=b(wrongCells>0);
    regCells=(frMean_HC(onOff>0,i)-frMean_spont(onOff>0,i))>3;
    c=frMean_spont(onOff>0,i);
    c=c(regCells>0);   
    
    ambiCells=(frMean_HC(onOff>0,i)-frMean_spont(onOff>0,i))>=-3&(frMean_HC(onOff>0,i)-frMean_spont(onOff>0,i))<=3;
    d=frMean_spont(onOff>0,i);
    d=d(ambiCells>0);  
    plot(c,'.b')     
    hold on
    plot(b,'.r')
    plot(d,'.g')
end




figure
nds='87654321'
col='rgb';
for i=1:8
    subplot(2,4,i)
    hold on
    k=frMean_spont(onOff>0&real_chirp'>0,i);

    a=mean(chirp(:,(i-1)*4+1:i*4,onOff>0&real_chirp'>0),2);
    a=reshape(a,23000,size(a,3));
    a=mean(a(500:2500,:));
    tmp=wf{i}(:,onOff>0&real_chirp'>0);
    tmp=mean(tmp(50:450,:));
    
    plot(k,'r')
    hold on
    plot(a,'b')
    plot(tmp,'g')

    
end
legend('HC','CH','QU')





figure
nds='87654321'
k=onOff>0;
for i=1:8
    subplot(4,2,i)
    hold on
    wrongCells=(frMean_HC(:,i)-frMean_spont(:,i))<-3;
    ambiCells=(frMean_HC(:,i)-frMean_spont(:,i))>=-3&(frMean_HC(:,i)-frMean_spont(:,i))<=3;
    regCells=(frMean_HC(:,i)-frMean_spont(:,i))>3;
    a=mean(chirp(:,(i-1)*4+1:i*4,k&real_chirp'>0&wrongCells'>0),2);
    a=reshape(a,23000,size(a,3));
%     a=a-repmat(mean(a(500:2500,:)),23000,1);    
    a=mean(a,2);
    [ind,lags]=xcorr(a(12000:19500),flicker(12000:19500),'coeff');
    plot(lags,ind,'r')
    line([0 0],[-1,1],'color','k')
    
    a=mean(chirp(:,(i-1)*4+1:i*4,k&real_chirp'>0&regCells'>0),2);
    a=reshape(a,23000,size(a,3));
%     a=a-repmat(mean(a(500:2500,:)),23000,1);   
    a=mean(a,2);
    [ind,lags]=xcorr(a(12000:19500),flicker(12000:19500),'coeff');
     plot(lags,ind)
    line([0 0],[-1,1],'color','k')
    title(['ND',nds(i),'   wrong cells ',int2str(sum(k&real_chirp'>0&wrongCells'>0)),'   regular cells ',int2str(sum(k&real_chirp'>0&regCells'>0)),'   ambi cells ',int2str(sum(k&real_chirp'>0&ambiCells'>0))])
end




plot(a)
hold on
plot(b-250:501:19523,a(b-250:501:19523),'*r')
plot(b:501:19523,a(b:501:19523),'*r')

for indd=1:length(maxpos)
    [maxx(indd) indmax(indd)]=max(a(maxpos(indd)-100:maxpos(indd)+100));
    indmax(indd)=indmax(indd)+maxpos(indd)-100-1;
    [minn(indd) indmin(indd)]=min(a(minneg(indd)-100:minneg(indd)+100));
    indmin(indd)=indmin(indd)+minneg(indd)-100-1;
end
plot(a)
hold on
plot(indmax,a(indmax),'*r')
plot(indmin,a(indmin),'*r')






figure
k=onOff>0;
for i=1:8
    subplot(2,4,i)
    wrongCells=(frMean_HC(k,i)-frMean_spont(k,i))<-3;
    regCells=(frMean_HC(k,i)-frMean_spont(k,i))>3;

    tmp=wf{i}(:,k);
%     tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(mean(tmp(:,wrongCells),2),'r');
    hold on
    plot(mean(tmp(:,regCells),2));
    hold on
    plot(mean(tmp(:,~regCells&~wrongCells),2),'g');
    title([int2str(sum(wrongCells)),' ',int2str(sum(regCells)),' ', int2str(sum(~regCells&~wrongCells))])
end



