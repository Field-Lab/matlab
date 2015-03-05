
clear

dates=cell(9,1);
dates{1}='20130227'
dates{2}='20130301'
dates{3}='20130301_1'
dates{4}='20130301_2'
dates{5}='20130302'
dates{6}='20130302_1'
dates{7}='20120329'
dates{8}='20120627'
dates{9}='20121023'

wfON=cell(8,1);
bfON=cell(8,1);
wfOFF=cell(8,1);
bfOFF=cell(8,1);
kmeanON=[];
kmeanOFF=[];

for datesCNT=1:9
    date=dates{datesCNT}
    path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
    load([path2save,'quick'],'white_flash','black_flash','names')
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date],'onOff','kmean','kstd','formula')
    kmeanON=[kmeanON kmean(formula,onOff>0)];
    kmeanOFF=[kmeanOFF kmean(formula,onOff<0)];
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



%% Plots



figure
subplot(2,1,1)
hold on
plot(kmeanON)
line([1 8],[0,0],'color','k','linewidth',3)
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
title('Firing Rate Difference: High contrast-low contrast ON CELLS')
xlabel('ND')

subplot(2,1,2)
hold on
plot(kmeanOFF)
line([1 8],[0,0],'color','k','linewidth',3)
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
title('Firing Rate Difference: High contrast-low contrast OFF CELLS')

hist([sum(kmeanON<0), 9])


i=2
a=find(kmeanON(i,:)<-5);

for m=1:size(a,2)
    subplot(5,3,m)
    for j=2:6
        plot(wfON{j}(:,a(m)))
        hold on
    end
end