
clear

dates=cell(9,1);
dates{1}='20130301'
dates{2}='20130301_1'
dates{3}='20130301_2'
dates{4}='20130302'
dates{5}='20130302_1'
dates{6}='20120329'
dates{7}='20121023'
dates{8}='20120627'
dates{9}='20130227'



zc_HC_ON=[];zc_LC_ON=[];
zc_HC_OFF=[];zc_LC_OFF=[];
wfON=cell(8,1);
bfON=cell(8,1);
wfOFF=cell(8,1);
bfOFF=cell(8,1);
frMean_HCON=[];
frMean_HCOFF=[];
frMean_LCON=[];
frMean_LCOFF=[];
frSTD_HCON=[];
frSTD_HCOFF=[];
frSTD_LCON=[];
frSTD_LCOFF=[];

frspont_ON=[];
frspont_OFF=[];
frSTDspont_ON=[];
frSTDspont_OFF=[];
nameON=[];
nameOFF=[];

for datesCNT=1:9
    date=dates{datesCNT}
    path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
    load([path2save,'quick'],'white_flash','black_flash','names')
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date],'onOff','zc_HC','zc_LC','frSTDspont','frMean_HC','frMean_LC','frSTD_HC','frSTD_LC','frspont','formula')
    nameON=[nameON; names(onOff>0)];
    nameOFF=[nameOFF; names(onOff<0)];
    zc_HC_ON=[zc_HC_ON; zc_HC(onOff>0,:)];
    zc_LC_ON=[zc_LC_ON; zc_LC(onOff>0,:)];
    zc_HC_OFF=[zc_HC_OFF; zc_HC(onOff<0,:)];
    zc_LC_OFF=[zc_LC_OFF; zc_LC(onOff<0,:)];
    frMean_HCON=[frMean_HCON frMean_HC(formula,onOff>0)];
    frMean_HCOFF=[frMean_HCOFF frMean_HC(formula,onOff<0)];
    frMean_LCON=[frMean_LCON frMean_LC(formula,onOff>0)];
    frMean_LCOFF=[frMean_LCOFF frMean_LC(formula,onOff<0)];
    frSTD_HCON=[frSTD_HCON frSTD_HC(formula,onOff>0)];
    frSTD_HCOFF=[frSTD_HCOFF frSTD_HC(formula,onOff<0)];
    frSTD_LCON=[frSTD_LCON frSTD_LC(formula,onOff>0)];
    frSTD_LCOFF=[frSTD_LCOFF frSTD_LC(formula,onOff<0)];
    frspont_ON=[frspont_ON frspont(formula,onOff>0)];
    frspont_OFF=[frspont_OFF frspont(formula,onOff<0)];
    frSTDspont_ON=[frSTDspont_ON frSTDspont(formula,onOff>0)];
    frSTDspont_OFF=[frSTDspont_OFF frSTDspont(formula,onOff<0)];
    size(black_flash,2)
    t=1;
    if size(black_flash,2)==72
        for i=8:9:72   
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
        for i=3:4:32
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
        for i=2:2:16
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
        for i=2:2:16
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
        for i=3:3:24
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
        for i=4:5:40
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
frMean_HCON=frMean_HCON';
frMean_HCOFF=frMean_HCOFF';
frMean_LCON=frMean_LCON';
frMean_LCOFF=frMean_LCOFF';
frspont_OFF=frspont_OFF';
frspont_ON=frspont_ON';
frSTDspont_OFF=frSTDspont_OFF';
frSTDspont_ON=frSTDspont_ON';
frSTD_HCON=frSTD_HCON';
frSTD_HCOFF=frSTD_HCOFF';
frSTD_LCON=frSTD_LCON';
frSTD_LCOFF=frSTD_LCOFF';

clear white_flash black_flash frSTDspont date datesCNT frspont formula frSTD_HC frSTD_LC frMean_HC frMean_LC names onOff path2save t tmp zc_HC zc_LC



%% Firing rate vs zero crossing


figure
nds='87654321'
for i=1:8
    subplot(4,2,i)
    plot(zc_HC_ON(:,i),frMean_HCON(:,i),'.r')
    hold on
    plot(zc_LC_ON(:,i),frMean_LCON(:,i),'.m')
%     plot(zc_HC_OFF(:,i),frMean_HCOFF(:,i),'.b')
%     plot(zc_LC_OFF(:,i),frMean_LCOFF(:,i),'.c')
    xlabel('zero crossing')
    ylabel('firing rate change')
    title(['ND',nds(i),'   red ON, blue OFF'])
    axis([70 250 0 80]);
    lims=get(gca);
    line([0,0],lims.YLim,'color','k')
    line(lims.XLim,[0 0],'color','k')
end


%% Firing rate at different contrasts
figure
nds='87654321'
for i=1:8
    subplot(3,3,i)
    plot(frMean_HCON(:,i)-frspont_ON(:,i),frMean_HCON(:,i)-frMean_LCON(:,i),'ro')
    hold on
    plot(frMean_HCOFF(:,i)-frspont_OFF(:,i),frMean_HCOFF(:,i)-frMean_LCOFF(:,i),'o')
    line([-20 40],[0,0],'color','k')
    line([0 0],[-20,40],'color','k')
    line([-20 40],[-20,40],'color','k')
    title(['Firing Rate changes at ND', nds(i)])
    xlabel('HC - spont')
    ylabel('HC - LC')
    axis([-20 40 -20 40])
end



figure
nds='87654321'
for i=1:8
    subplot(3,3,i)
    plot(frMean_HCON(:,i)-frspont_ON(:,i),frMean_LCON(:,i)-frspont_ON(:,i),'ro')
    hold on
    plot(frMean_HCOFF(:,i)-frspont_OFF(:,i),frMean_LCOFF(:,i)-frspont_OFF(:,i),'o')
    line([-20 40],[0,0],'color','k')
    line([0 0],[-20,40],'color','k')
    line([-20 40],[-20,40],'color','k')
    title(['Firing Rate changes at ND', nds(i)])
    xlabel('HC - spont')
    ylabel('LC - spont')
    axis([-20 40 -20 40])
end



figure
nds='87654321'
for i=1:8
    subplot(3,3,i)
    plot(frSTD_HCON(:,i)-frSTDspont_ON(:,i),frSTD_HCON(:,i)-frSTD_LCON(:,i),'ro')
    hold on
    plot(frSTD_HCOFF(:,i)-frSTDspont_OFF(:,i),frSTD_HCOFF(:,i)-frSTD_LCOFF(:,i),'o')
    line([-20 40],[0,0],'color','k')
    line([0 0],[-20,40],'color','k')
    line([-20 40],[-20,40],'color','k')
    title(['Firing Rate changes at ND', nds(i)])
    xlabel('HC - spont')
    ylabel('HC - LC')
    axis([-20 40 -20 40])
end


figure
nds='87654321'
for i=1:8
    subplot(3,3,i)
    plot(frSTD_LCON(:,i)-frSTDspont_ON(:,i),frSTD_HCON(:,i)-frSTDspont_ON(:,i),'ro')
    hold on
    plot(frSTD_LCOFF(:,i)-frSTDspont_OFF(:,i),frSTD_HCOFF(:,i)-frSTDspont_OFF(:,i),'o')
    line([-20 40],[0,0],'color','k')
    line([0 0],[-20,40],'color','k')
    line([-20 40],[-20,40],'color','k')
    title(['Firing Rate changes: variance at ND', nds(i)])
    xlabel('LC - spont')
    ylabel('HC - spont')
    axis([-10 30 -10 30])
end



k=sum(a>=-2&a<=2)/size(a,1)*100;
k1=sum(a>2&a<=10)/size(a,1)*100;
k2=sum(a<-2&a>-10)/size(a,1)*100;
figure
subplot(2,1,1)
bar([k1; k; k2]')
axis([0 9 0 50])
legend({'speed up','no change', 'slow down'})
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('% of total')
title('ON cells: changes at low contrast')
a=zc_HC_OFF-zc_LC_OFF;
k=sum(a>=-2&a<=2)/size(a,1)*100;
k1=sum(a>2&a<=10)/size(a,1)*100;
k2=sum(a<-2&a>-10)/size(a,1)*100;
subplot(2,1,2)
bar([k1; k; k2]')
axis([0 9 0 50])
legend({'speed up','no change', 'slow down'})
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('% of total')
title('OFF cells: changes at low contrast')


%% Firing rate vs quick

i=3
wrongCells=(frMean_HCON(:,i)-frspont_ON(:,i))<-2;
wrongOffCells=(frMean_HCOFF(:,i)-frspont_OFF(:,i))<-2;

figure
for i=1:8
    subplot(2,4,i)
    a=wfON{i}(:,~wrongCells);
    a=a-repmat(mean(a(1:490,:)),4500,1);
    hold on
    plot(mean(a,2),'b','linewidth',2)
    a=bfON{i}(:,~wrongCells);
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(4701:9200,mean(a,2),'b','linewidth',2)
    
    a=wfON{i}(:,wrongCells);
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(mean(a,2),'r','linewidth',2)
    a=bfON{i}(:,wrongCells);
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(4701:9200,mean(a,2),'r','linewidth',2)
    
    axis([0 9200 -30 50])
end

figure
for i=2:7
    subplot(2,3,i-1)
    a=wfON{i}(:,~wrongCells);
    a=a-repmat(mean(a(1:490,:)),4500,1);
    plot(a)
    hold on
    plot(mean(a,2),'k','linewidth',2)
end


i=3
b=find(wrongCells);
for j=1:sum(wrongCells)
    subplot(6,6,j)
    a=wfON{i}(:,b(j));
    a=a-mean(a(1:490,:));
    plot(a)
    title([nameON{b(j)}, '  ', num2str(frMean_HCON(b(j),i)-frMean_LCON(b(j),i))],'Interpreter','none')
end


wrongCells=(frMean_HCON(:,2:6)-frMean_LCON(:,2:6))<-2;
figure
hist([sum(wrongCells,2);9])

%% Linear Filters


clear

dates=cell(9,1);
dates{1}='20130301'
dates{2}='20130301_1'
dates{3}='20130301_2'
dates{4}='20130302'
dates{5}='20130302_1'
dates{6}='20120329'
dates{7}='20121023'
dates{8}='20120627'
dates{9}='20130227'
HC=[];
LC=[];
onOff_all=[];
for datesCNT=1:9
    date=dates{datesCNT}
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date],'onOff','LinearFilter')
    onOff_all=[onOff_all onOff];
    t=1;
    tt=size(HC,3)+1*~isempty(HC):size(LinearFilter,4)+size(HC,3)-1*isempty(HC);
    if size(LinearFilter,3)==16
        for i=1:2:16
            HC(:,t,tt)=reshape(mean(mean(LinearFilter(:,1:2:end,i:i+1,:),2),3),500,size(LinearFilter,4));
            LC(:,t,tt)=reshape(mean(mean(LinearFilter(:,2:2:end,i:i+1,:),2),3),500,size(LinearFilter,4));
            t=t+1;
        end
    elseif size(LinearFilter,3)==24
        for i=1:3:24
            HC(:,t,tt)=reshape(mean(mean(LinearFilter(:,1:2:end,i:i+2,:),2),3),500,size(LinearFilter,4));
            LC(:,t,tt)=reshape(mean(mean(LinearFilter(:,2:2:end,i:i+2,:),2),3),500,size(LinearFilter,4));
            t=t+1;
        end
    elseif size(LinearFilter,3)==46
        for i=1:4:32
            HC(:,t,tt)=reshape(mean(mean(LinearFilter(:,1:2:end,i:i+3,:),2),3),500,size(LinearFilter,4));
            LC(:,t,tt)=reshape(mean(mean(LinearFilter(:,2:2:end,i:i+3,:),2),3),500,size(LinearFilter,4));
            t=t+1;
        end
    elseif size(LinearFilter,3)==342
        for i=1:24:24*8
            HC(:,t,tt)=reshape(mean(LinearFilter(:,1,i:2:i+23,:),3),500,size(LinearFilter,4));
            LC(:,t,tt)=reshape(mean(LinearFilter(:,1,(i+1):2:i+23,:),3),500,size(LinearFilter,4));
            t=t+1;
        end
    elseif size(LinearFilter,3)==49
        for i=1:6:6*8
            HC(:,t,tt)=reshape(mean(mean(LinearFilter(:,1:2:end,i:i+5,:),2),3),500,size(LinearFilter,4));
            LC(:,t,tt)=reshape(mean(mean(LinearFilter(:,2:2:end,i:i+5,:),2),3),500,size(LinearFilter,4));
            t=t+1;
        end
    elseif size(LinearFilter,3)==36
        for i=5:4:36
            HC(:,t,tt)=reshape(mean(mean(LinearFilter(:,1:2:end,i:i+3,:),2),3),500,size(LinearFilter,4));
            LC(:,t,tt)=reshape(mean(mean(LinearFilter(:,2:2:end,i:i+3,:),2),3),500,size(LinearFilter,4));
            t=t+1;
        end
    end    

end
onOff=onOff_all;
clear onOff_all LinearFilter date dates datesCNT t tt
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/LinearFilter9exp','onOff','HC','LC')

%% Plots
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/LinearFilter9exp','onOff','HC','LC')
HC_ON=HC(:,:,onOff>0);
HC_OFF=HC(:,:,onOff<0);
LC_ON=LC(:,:,onOff>0);
LC_OFF=LC(:,:,onOff<0);

frMean_HCON=frMean_HC(onOff>0,:);
frspont_ON=frMean_spont(onOff>0,:);
i=5
wrongCells=(frMean_HCON(:,i)-frspont_ON(:,i))<-2;

figure
for i=1:8
    subplot(2,4,i)
    t=reshape(LC_ON(:,i,~wrongCells),500,sum(~wrongCells));
%     t=t-repmat(mean(t),size(t,1),1);
    t=t./repmat(sum(abs(t)),size(t,1),1);
    plot(nanmean(t,2),'r')
%     plot(t,'r')
    hold on
    t=reshape(LC_ON(:,i,wrongCells),500,sum(wrongCells));
%     t=t-repmat(mean(t),size(t,1),1);
    t=t./repmat(sum(abs(t)),size(t,1),1);
%     plot(t,'b')
    plot(nanmean(t,2))
    title('Low Contrast')
    legend({'st cells','wrong cells'})
    line([0 500],[0,0],'color','k')

end

figure

k=wrongCells;
for i=1:8
    subplot(2,4,i)
    t=reshape(HC(:,i,onOff<0),500,sum(onOff<0));
%     t=t-repmat(mean(t),size(t,1),1);
%     t=t./repmat(sum(abs(t)),size(t,1),1);
    plot(nanmean(t,2))
    hold on
    t=reshape(HC_ON(:,i,k),500,sum(k));
%     t=t-repmat(mean(t),size(t,1),1);
%     t=t./repmat(sum(abs(t)),size(t,1),1);
    plot(-nanmean(t,2),'r')
    
    t=reshape(HC_ON(:,i,~k),500,sum(~k));
    %     t=t-repmat(mean(t),size(t,1),1);
    %     t=t./repmat(sum(abs(t)),size(t,1),1);
    plot(-nanmean(t,2),'g')
    legend({'OFF','ON wrong','ON right'})
end





figure
nds='87654321'
for i=1:8
    subplot(2,4,i)
    k=(frMean_HCON(:,i)-frspont_ON(:,i))<-0.2*mean(frspont_ON(:,i));
    t=reshape(LC_ON(:,i,k),500,sum(k));
%     t=t-repmat(mean(t),size(t,1),1);
    t=t./repmat(sum(abs(t)),size(t,1),1);
%     plot(t,'b')
    plot(nanmean(t,2))
    hold on
    t=reshape(HC_ON(:,i,k),500,sum(k));
%     t=t-repmat(mean(t),size(t,1),1);
    t=t./repmat(sum(abs(t)),size(t,1),1);
%     plot(t,'b')
    plot(nanmean(t,2),'r')
    legend({'ON wrong low','ON wrong high'})
    line([0,500],[0,0],'color','k')
    title([nds(i),'  ',int2str(sum(k))])
end



figure
nds='87654321'
for i=1:8
    subplot(2,4,i)
    k=(frMean_HCON(:,i)-frspont_ON(:,i))>0.2*mean(frspont_ON(:,i));
    t=reshape(LC_ON(:,i,k),500,sum(k));
%     t=t-repmat(mean(t),size(t,1),1);
    t=t./repmat(sum(abs(t)),size(t,1),1);
%     plot(t,'b')
    plot(nanmean(t,2))
    hold on
    t=reshape(HC_ON(:,i,k),500,sum(k));
%     t=t-repmat(mean(t),size(t,1),1);
    t=t./repmat(sum(abs(t)),size(t,1),1);
%     plot(t,'b')
    plot(nanmean(t,2),'r')
    legend({'ON reg low','ON reg high'})
    line([0,500],[0,0],'color','k')
    title([nds(i),'  ',int2str(sum(k))])
end







figure
for i=1:8
    subplot(2,4,i)
    t=reshape(LC_OFF(:,i,:),500,sum(onOff<0));
%     t=t-repmat(mean(t),size(t,1),1);
    t=t./repmat(sum(abs(t)),size(t,1),1);
%     plot(t,'b')
    plot(nanmean(t,2))
    hold on
    t=reshape(HC_OFF(:,i,:),500,sum(onOff<0));
%     t=t-repmat(mean(t),size(t,1),1);
    t=t./repmat(sum(abs(t)),size(t,1),1);
%     plot(t,'b')
    plot(nanmean(t,2),'r')
    line([0,500],[0,0],'color','k')

end