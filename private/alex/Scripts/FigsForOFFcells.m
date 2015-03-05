%% Accumulate

clear

dates=cell(12,1);
dates{1}='20130220'
dates{2}='20130220_1'
dates{3}='20130224'
dates{4}='20130225'
dates{5}='20130226'
dates{6}='20130227'
dates{7}='20130301'
dates{8}='20130301_1'
dates{9}='20130301_2'
dates{10}='20130302'
dates{11}='20130302_1'
dates{12}='20120329'


wfON=cell(8,1);
bfON=cell(8,1);
wfOFF=cell(8,1);
bfOFF=cell(8,1);
nameON=[];
nameOFF=[];

for datesCNT=1:12
    date=dates{datesCNT}
    path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
    load([path2save,'quick'],'white_flash','black_flash','names')
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date],'onOff')
    size(black_flash,2)

    nameON=[nameON; names(onOff>0)];
    nameOFF=[nameOFF; names(onOff<0)];
    t=1;
    if size(black_flash,2)==72
        for i=1:9:72   
            tmp=mean(white_flash(:,i:i+1,onOff>0),2);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=mean(black_flash(:,i:i+1,onOff>0),2);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=mean(white_flash(:,i:i+1,onOff<0),2);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=mean(black_flash(:,i:i+1,onOff<0),2);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
    elseif size(black_flash,2)==32
        for i=1:4:32
            tmp=mean(white_flash(:,i:i+1,onOff>0),2);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=mean(black_flash(:,i:i+1,onOff>0),2);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=mean(white_flash(:,i:i+1,onOff<0),2);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=mean(black_flash(:,i:i+1,onOff<0),2);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
    elseif size(black_flash,2)==16
        for i=1:2:16
            tmp=white_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=black_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=white_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=black_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
    elseif size(black_flash,2)==17
        for i=1:2:16
            tmp=white_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=black_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=white_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=black_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
    elseif size(black_flash,2)==24
        for i=1:3:24
            tmp=white_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=black_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=white_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=black_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
    elseif size(black_flash,2)==56
        for i=1:5:40
            tmp=mean(white_flash(:,i:i+1,onOff>0),2);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=mean(black_flash(:,i:i+1,onOff>0),2);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=mean(white_flash(:,i:i+1,onOff<0),2);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=mean(black_flash(:,i:i+1,onOff<0),2);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
    elseif size(black_flash,2)==8
        for i=1:8
            tmp=white_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=black_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=white_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=black_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
    elseif size(black_flash,2)==7
        t=2;
        for i=1:7
            tmp=white_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=black_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=white_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=black_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
        wfON{1}=[wfON{1} zeros(4500,sum(onOff>0))];
        wfOFF{1}=[wfOFF{1} zeros(4500,sum(onOff<0))];
        bfON{1}=[bfON{1} zeros(4500,sum(onOff>0))];
        bfOFF{1}=[bfOFF{1} zeros(4500,sum(onOff<0))];
    elseif size(black_flash,2)==14
        for i=1:8
            tmp=white_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            wfON{t}=[wfON{t} tmp];
            
            tmp=black_flash(:,i,onOff>0);
            tmp=reshape(tmp,4500,sum(onOff>0));
            bfON{t}=[bfON{t} tmp];
            
            tmp=white_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            wfOFF{t}=[wfOFF{t} tmp];
            
            tmp=black_flash(:,i,onOff<0);
            tmp=reshape(tmp,4500,sum(onOff<0));
            bfOFF{t}=[bfOFF{t} tmp];
            t=t+1;
        end
    end    

end

clear black_flash white_flash date dates dateCNT tmp path2save names onOff



figure
nds='87654321'
for i=1:8
    subplot(4,2,i)
    tmp=wfOFF{i};
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    hold off
    plot(mean(tmp,2),'color','k','linewidth',2)
    %indices
%     [a,b]=max(tmp(ind_start1-100:ind_start1+100,:));
%     ampl1(i,1)=mean(a);
%     ampl1(i,2)=std(a)/sqrt(length(a));
% 
%     lat1(i,1)=mean(b);
%     lat1(i,2)=std(b)/sqrt(length(b));
%     
%     
    hold on
    tmp=bfOFF{i};
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(4701:9200,mean(tmp,2),'color','k','linewidth',2)

    %indices
%     [a,b]=max(tmp(ind_start2-100:ind_start2+100,:));
%     ampl2(i,1)=mean(a);
%     ampl2(i,2)=std(a)/sqrt(length(a));
% 
%     lat2(i,1)=mean(b);
%     lat2(i,2)=std(b)/sqrt(length(b));    
%     
    
%     
%     axis([1 9200 -5 65])
%     line([1 4500],[0,0],'color','k')
%     line([4700 9200],[0,0],'color','k')
%     
%     line([1 500],[55,55],'color','k')
%     line([500 2500],[60,60],'color','k')
%     line([2500 4500],[55,55],'color','k')
%     line([500 500],[55,60],'color','k')
%     line([2500 2500],[55,60],'color','k')
%     
%     line([1 500]+4700,[55,55],'color','k')
%     line([500 2500]+4700,[50,50],'color','k')
%     line([2500 4500]+4700,[55,55],'color','k')
%     line([500 500]+4700,[50,55],'color','k')
%     line([2500 2500]+4700,[50,55],'color','k')


% 
    axis([1 9200 -20 50])
    line([1 4500],[0,0],'color','k')
    line([4700 9200],[0,0],'color','k')
    
    line([1 500],[40,40],'color','k')
    line([500 2500],[45,45],'color','k')
    line([2500 4500],[40,40],'color','k')
    line([500 500],[40,45],'color','k')
    line([2500 2500],[40,45],'color','k')
    
    line([1 500]+4700,[40,40],'color','k')
    line([500 2500]+4700,[35,35],'color','k')
    line([2500 4500]+4700,[40,40],'color','k')
    line([500 500]+4700,[35,40],'color','k')
    line([2500 2500]+4700,[35,40],'color','k')
    title(['ND',nds(i)])
end



figure
nds='87654321'
for i=1:8
    subplot(4,2,i)
    tmp=wf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);

    plot(mean(tmp,2),'color','r','linewidth',2)
 
    hold on
    tmp=bf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(4701:9200,mean(tmp,2),'color','r','linewidth',2)



    axis([1 9200 -20 50])
    line([1 4500],[0,0],'color','k')
    line([4700 9200],[0,0],'color','k')
    
    line([1 500],[40,40],'color','k')
    line([500 2500],[45,45],'color','k')
    line([2500 4500],[40,40],'color','k')
    line([500 500],[40,45],'color','k')
    line([2500 2500],[40,45],'color','k')
    
    line([1 500]+4700,[40,40],'color','k')
    line([500 2500]+4700,[35,35],'color','k')
    line([2500 4500]+4700,[40,40],'color','k')
    line([500 500]+4700,[35,40],'color','k')
    line([2500 2500]+4700,[35,40],'color','k')
    title(['ND',nds(i)])
end




nds='87654321'
for i=1:8
    subplot(4,2,i)
    tmp=bf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:200
        [k(i,j),m]=max(tmp(500:1000,j));
        b=find(tmp(500+m:2500,j)<a(j),1);
        if ~isempty(b)
            returns(j,i)=b+m;
        else
            returns(j,i)=2000;
        end
    end

    tmp=wf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:200
        [kWF(i,j),m]=max(tmp(2500:3000,j));
        b=find(tmp(2500+m:4500,j)<a(j),1);
        if ~isempty(b)
            returnsWF(j,i)=b+m;
        else
            returnsWF(j,i)=2000;
        end
    end
    medWF(i)=mean(returnsWF(returnsWF(:,i)<1000,i));
    medBF(i)=mean(returns(returns(:,i)<1000,i));

end


nds='87654321'
for i=1:8
    subplot(4,2,i)
    tmp=bf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:200
        [~,m]=max(tmp(500:1000,j));
        b=find(tmp(500+m:2500,j)<a(j),1);
        if ~isempty(b)
            returns(j,i)=b+m;
        else
            returns(j,i)=2000;
        end
    end

    tmp=wf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:200
        [~,m]=max(tmp(2500:3000,j));
        b=find(tmp(2500+m:4500,j)<a(j),1);
        if ~isempty(b)
            returnsWF(j,i)=b+m;
        else
            returnsWF(j,i)=2000;
        end
    end

    hist([returns(:,i)-returnsWF(:,i); -2000; 2000],80)
    axis([-2000 2000 0 30])
%     plot(returns(:,i),returnsWF(:,i),'o')
%     xlabel('Black latency')
%     ylabel('White latency')
    title(['ND',nds(i),' black - white'])
end




figure
errorbar(mean(returns),std(returns)/sqrt(200))
hold on
errorbar(mean(returnsWF),std(returnsWF)/sqrt(200),'r')
for i=1:8
    [h,p(i)]=ttest(returns(:,i),returnsWF(:,i));
    if h
        plot(i, mean(returns(:,i))+100,'*k')
    end
end
xlabel('ND')
ylabel('HZ')
title('Mean +/-st error of latency to reach mean+std of the spontaneous firing rate in OFF response')
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
legend('Negative step','Positive step','location','southeast')




figure
subplot(1,2,1)
bar([sum(returns>1999); sum(returnsWF>1999)]'/2)
xlabel('ND')
ylabel('HZ')
title('% of cells not returning to the baseline in 2s')
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
legend('Negative step','Positive step')
subplot(1,2,2)
bar([sum(returns<300); sum(returnsWF<300)]'/2)
xlabel('ND')
ylabel('% of all OFF cells')
title('% of cells returning to the baseline in 600ms')
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
legend('Negative step','Positive step')


figure
hist([returns(:); returnsWF(:)],50)





figure
tt=[3,5]
subplot(1,2,1)
hist([max(returns(:,tt)')-min(returns(:,tt)')],4);
subplot(1,2,2)
hist([max(returnsWF(:,tt)')-min(returnsWF(:,tt)')],4);



for i=1:8
    tmp=bf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);    
    tmp1=wf{i}(:,onOff<0);
    tmp1=tmp1-repmat(mean(tmp1(50:450,:)),4500,1);
    for j=1:200
        cors(i,j)=corr(tmp(501:2500,j),tmp1(2501:4500,j));
        corsWF(i,j)=corr(tmp(2501:4500,j),tmp1(501:2500,j));
    end
end
figure
bar([sum(cors'<-0.2)/2; sum(corsWF'<-0.2)/2]')
xlabel('ND')
ylabel('% of all OFF cells')
title('% of cell with correlation of OFF / ON responses <-0.2')
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
axis([0 9 0 14])
legend('OFF','ON')

figure
i=4;
p=find(corsWF(i,:)<-0.2);
for j=1:sum(corsWF(i,:)<-0.2)
    subplot(6,5,j)
    tmp=bf{i}(:,onOff<0); 
    tmp1=wf{i}(:,onOff<0);
    plot(tmp(:,p(j)))
    hold on
    plot(4701:9200,tmp1(:,p(j)),'r')
end





