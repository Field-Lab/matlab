clear
addpath(genpath('/Users/alexth/Desktop/Scripts'))
path2load=['/Users/alexth/Desktop/old_stuff/ONresp/'];
path2save=['/Users/alexth/Desktop/old_stuff/ONresp/allDates/'];
allDates=dir([path2load,'2013*']);

for i=1:length(allDates)
    a=dir([path2load,allDates(i).name,'/quick_plot/*.png']);
    for j=1:length(a)
        copyfile([path2load,allDates(i).name,'/quick_plot/',a(j).name],[path2save,a(j).name])
    end
end




%% Cumulative plot

clear

dates{1}='2013-09-03';
dates{2}='20130723a';
dates{3}='20130723b';
dates{4}='20130906';

for ddd=1:4
    
    date=dates{ddd};
    path2load=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
    load([path2load,'quick'])
    load([path2load,'LF_HC'])
    
    path2save=['/Users/alexth/Desktop/old_stuff/ONresp/cplot/'];
    if ~exist(path2save,'dir')
        mkdir(path2save);
    end
    
    
    figure
    set(gcf,'position',[4          87        1243         849])
    
    if strcmp(date,'20130723a')
        startFrom=[3 7 15 24 28 36 47 51 59 71];
    elseif strcmp(date,'2013-09-03')
        startFrom=[3 7 15 27 31 39 51 55 63 77];
    else
        startFrom=[3 7 15 27 31 39 51 55 63 75]; % number of trial to start averaging. First is for ND7, then 3 for ND6 (control, APB, wash), 3 fopr ND5 and 3 for ND4. Take 4 trials for each average from the middle of the run.
    end
    
    col='brg';tit='7654';sbpl=[3 5 7];
    
    for j=1:size(black_flash,3)
        
        % plot quick
        subplot(4,2,1)
        hold off
        a=[mean(black_flash(:,startFrom(1):startFrom(1)+3,j),2); zeros(200,1); mean(white_flash(:,startFrom(1):startFrom(1)+3,j),2)];
        plot(a,'b','linewidth',2)
        hold on
        %plot black flash timing
        line([0 500], [max(a)*1.2, max(a)*1.2],'color','k')
        line([500 500], [max(a)*1.1, max(a)*1.2],'color','k')
        line([500 2500], [max(a)*1.1, max(a)*1.1],'color','k')
        line([2500 2500], [max(a)*1.1, max(a)*1.2],'color','k')
        line([2500 4500], [max(a)*1.2, max(a)*1.2],'color','k')
        %plot white flash timing
        line([4700 5300], [max(a)*1.2, max(a)*1.2],'color','k')
        line([5300 5300], [max(a)*1.2, max(a)*1.3],'color','k')
        line([5300 7300], [max(a)*1.3, max(a)*1.3],'color','k')
        line([7300 7300], [max(a)*1.2, max(a)*1.3],'color','k')
        line([7300 9300], [max(a)*1.2, max(a)*1.2],'color','k')
        %fill gap
        for i=1:15
            line([4500 4700], [max(a)*1.4/15*(i-1),max(a)*1.4/15*(i)],'color','k')
        end
        title('ND7')
        axis([0 9200 0 max(max(a)*1.35,1)])
        
        k=1;
        for cnt=[2,5,8]
            subplot(4,2,sbpl(k))
            hold off
            a=[mean(black_flash(:,startFrom(cnt):startFrom(cnt)+3,j),2); zeros(200,1); mean(white_flash(:,startFrom(cnt):startFrom(cnt)+3,j),2)];
            plot(a,col(1),'linewidth',2)
            tmp=max(a);
            hold on
            a=[mean(black_flash(:,startFrom(cnt+1):startFrom(cnt+1)+3,j),2); zeros(200,1); mean(white_flash(:,startFrom(cnt+1):startFrom(cnt+1)+3,j),2)];
            plot(a,col(2),'linewidth',2)
            tmp=max(a,tmp);
            a=[mean(black_flash(:,startFrom(cnt+2):startFrom(cnt+2)+3,j),2); zeros(200,1); mean(white_flash(:,startFrom(cnt+2):startFrom(cnt+2)+3,j),2)];
            plot(a,col(3),'linewidth',2)
            a=max(a,tmp);
            
            %plot black flash timing
            line([0 500], [max(a)*1.2, max(a)*1.2],'color','k')
            line([500 500], [max(a)*1.1, max(a)*1.2],'color','k')
            line([500 2500], [max(a)*1.1, max(a)*1.1],'color','k')
            line([2500 2500], [max(a)*1.1, max(a)*1.2],'color','k')
            line([2500 4500], [max(a)*1.2, max(a)*1.2],'color','k')
            %plot white flash timing
            line([4700 5300], [max(a)*1.2, max(a)*1.2],'color','k')
            line([5300 5300], [max(a)*1.2, max(a)*1.3],'color','k')
            line([5300 7300], [max(a)*1.3, max(a)*1.3],'color','k')
            line([7300 7300], [max(a)*1.2, max(a)*1.3],'color','k')
            line([7300 9300], [max(a)*1.2, max(a)*1.2],'color','k')
            %fill gap
            for i=1:15
                line([4500 4700], [max(a)*1.4/15*(i-1),max(a)*1.4/15*(i)],'color','k')
            end
            title(['ND',tit(k+1)])
            axis([0 9200 0 max(max(a)*1.35,1)])
            k=k+1;
        end
        
        
        
        % plot LF
        
        subplot(4,2,2)
        plot(mean(LinearFilter(:,startFrom(1):startFrom(1)+3,j),2),'b','linewidth',2)
        line([0 500],[0,0],'color','k')
        title('ND7')
        
        
        
        k=1;
        for cnt=[2,5,8]
            subplot(4,2,sbpl(k)+1)
            
            hold off
            plot(mean(LinearFilter(:,startFrom(cnt):startFrom(cnt)+3,j),2),col(1),'linewidth',2)
            hold on
            plot(mean(LinearFilter(:,startFrom(cnt+1):startFrom(cnt+1)+3,j),2),col(2),'linewidth',2)
            plot(mean(LinearFilter(:,startFrom(cnt+2):startFrom(cnt+2)+3,j),2),col(3),'linewidth',2)
            if cnt==8
                legend('control','APB','wash')
            end
            title(['ND',tit(k+1)])
            line([0 500],[0,0],'color','k')
            k=k+1;
        end
        
        
        subplot('Position',[0.5 0.95 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title([names{j},' Mean of 4 trials in the end of the condition block'],'FontSize',12,'FontWeight','bold','Interpreter','None')
        saveas(gcf,[path2save,names{j},'.png'])
    end
    close all
end






%% Cumulative data - averaged (checked for discarded units) by condition

clear

dates{1}='2013-09-03';
dates{2}='20130723a';
dates{3}='20130723b';
dates{4}='20130906';

datarun=struct('date',{},'name',{},'white',{},'black',{},'lf',{},'chirp',{});

walk=1;
for ddd=1:4
    
    date=dates{ddd};
    path2load=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
    load([path2load,'quick'])
    load([path2load,'LF_HC'])
    load([path2load,'chirp'])
    
    path2check=['/Users/alexth/Desktop/old_stuff/ONresp/cplot/'];
    names2check=dir(path2check);
    
    if strcmp(date,'20130723a')
        startFrom=[3 7 15 24 28 36 47 51 59 71];
    elseif strcmp(date,'2013-09-03')
        startFrom=[3 7 15 27 31 39 51 55 63 77];
    else
        startFrom=[3 7 15 27 31 39 51 55 63 75]; % number of trial to start averaging. First is for ND7, then 3 for ND6 (control, APB, wash), 3 fopr ND5 and 3 for ND4. Take 4 trials for each average from the middle of the run.
    end
    
    for i=1:size(black_flash,3)
        
        j=1;
        while ~strcmp(names{i},names2check(j).name(1:end-4)) && j<length(names2check)
            j=j+1;
        end        
        if strcmp(names{i},names2check(j).name(1:end-4))
            datarun(walk).date=date;
            datarun(walk).name=names{i};
            
            datarun(walk).black.nd7=mean(black_flash(:,startFrom(1):startFrom(1)+3,i),2);
            datarun(walk).black.control=[mean(black_flash(:,startFrom(2):startFrom(2)+3,i),2) mean(black_flash(:,startFrom(5):startFrom(5)+3,i),2) mean(black_flash(:,startFrom(8):startFrom(8)+3,i),2)];
            datarun(walk).black.apb=[mean(black_flash(:,startFrom(3):startFrom(3)+3,i),2) mean(black_flash(:,startFrom(6):startFrom(6)+3,i),2) mean(black_flash(:,startFrom(9):startFrom(9)+3,i),2)];
            datarun(walk).black.wash=[mean(black_flash(:,startFrom(4):startFrom(4)+3,i),2) mean(black_flash(:,startFrom(7):startFrom(7)+3,i),2) mean(black_flash(:,startFrom(10):startFrom(10)+3,i),2)];
            datarun(walk).white.nd7=mean(white_flash(:,startFrom(1):startFrom(1)+3,i),2);
            datarun(walk).white.control=[mean(white_flash(:,startFrom(2):startFrom(2)+3,i),2) mean(white_flash(:,startFrom(5):startFrom(5)+3,i),2) mean(white_flash(:,startFrom(8):startFrom(8)+3,i),2)];
            datarun(walk).white.apb=[mean(white_flash(:,startFrom(3):startFrom(3)+3,i),2) mean(white_flash(:,startFrom(6):startFrom(6)+3,i),2) mean(white_flash(:,startFrom(9):startFrom(9)+3,i),2)];
            datarun(walk).white.wash=[mean(white_flash(:,startFrom(4):startFrom(4)+3,i),2) mean(white_flash(:,startFrom(7):startFrom(7)+3,i),2) mean(white_flash(:,startFrom(10):startFrom(10)+3,i),2)];

            datarun(walk).chirp.nd7=mean(chirp(:,startFrom(1):startFrom(1)+3,i),2);
            datarun(walk).chirp.control=[mean(chirp(:,startFrom(2):startFrom(2)+3,i),2) mean(chirp(:,startFrom(5):startFrom(5)+3,i),2) mean(chirp(:,startFrom(8):startFrom(8)+3,i),2)];
            datarun(walk).chirp.apb=[mean(chirp(:,startFrom(3):startFrom(3)+3,i),2) mean(chirp(:,startFrom(6):startFrom(6)+3,i),2) mean(chirp(:,startFrom(9):startFrom(9)+3,i),2)];
            datarun(walk).chirp.wash=[mean(chirp(:,startFrom(4):startFrom(4)+3,i),2) mean(chirp(:,startFrom(7):startFrom(7)+3,i),2) mean(chirp(:,startFrom(10):startFrom(10)+3,i),2)];

            
            datarun(walk).lf.nd7=mean(LinearFilter(:,startFrom(1):startFrom(1)+3,i),2);
            datarun(walk).lf.control=[mean(LinearFilter(:,startFrom(2):startFrom(2)+3,i),2) mean(LinearFilter(:,startFrom(5):startFrom(5)+3,i),2) mean(LinearFilter(:,startFrom(8):startFrom(8)+3,i),2)];
            datarun(walk).lf.apb=[mean(LinearFilter(:,startFrom(3):startFrom(3)+3,i),2) mean(LinearFilter(:,startFrom(6):startFrom(6)+3,i),2) mean(LinearFilter(:,startFrom(9):startFrom(9)+3,i),2)];
            datarun(walk).lf.wash=[mean(LinearFilter(:,startFrom(4):startFrom(4)+3,i),2) mean(LinearFilter(:,startFrom(7):startFrom(7)+3,i),2) mean(LinearFilter(:,startFrom(10):startFrom(10)+3,i),2)];
           
            
            walk=walk+1;
        end
    end    
       
end

save('/Users/alexth/Desktop/old_stuff/ONresp/datarun.mat','datarun')



%% Cumulative data - averaged (checked for discarded units) by nd

clear

dates{1}='2013-09-03';
dates{2}='20130723a';
dates{3}='20130723b';
dates{4}='20130906';

datarun=struct('date',{},'name',{},'white',{},'black',{},'lf',{},'chirp',{});

walk=1;
for ddd=1:4
    
    date=dates{ddd};
    path2load=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
    load([path2load,'quick'])
    load([path2load,'LF_HC'])
    load([path2load,'chirp'])
    
    path2check=['/Users/alexth/Desktop/old_stuff/ONresp/cplot/'];
    names2check=dir(path2check);
    
    if strcmp(date,'20130723a')
        startFrom=[3 7 15 24 28 36 47 51 59 71];
    elseif strcmp(date,'2013-09-03')
        startFrom=[3 7 15 27 31 39 51 55 63 77];
    else
        startFrom=[3 7 15 27 31 39 51 55 63 75]; % number of trial to start averaging. First is for ND7, then 3 for ND6 (control, APB, wash), 3 fopr ND5 and 3 for ND4. Take 4 trials for each average from the middle of the run.
    end
    
    for i=1:size(black_flash,3)
        
        j=1;
        while ~strcmp(names{i},names2check(j).name(1:end-4)) && j<length(names2check)
            j=j+1;
        end        
        if strcmp(names{i},names2check(j).name(1:end-4))
            datarun(walk).date=date;
            datarun(walk).name=names{i};
            
            datarun(walk).black.nd7=mean(black_flash(:,startFrom(1):startFrom(1)+3,i),2);
            datarun(walk).black.nd6=[mean(black_flash(:,startFrom(2):startFrom(2)+3,i),2) mean(black_flash(:,startFrom(3):startFrom(3)+3,i),2) mean(black_flash(:,startFrom(4):startFrom(4)+3,i),2)];
            datarun(walk).black.nd5=[mean(black_flash(:,startFrom(5):startFrom(5)+3,i),2) mean(black_flash(:,startFrom(6):startFrom(6)+3,i),2) mean(black_flash(:,startFrom(7):startFrom(7)+3,i),2)];
            datarun(walk).black.nd4=[mean(black_flash(:,startFrom(8):startFrom(8)+3,i),2) mean(black_flash(:,startFrom(9):startFrom(9)+3,i),2) mean(black_flash(:,startFrom(10):startFrom(10)+3,i),2)];
            datarun(walk).white.nd7=mean(white_flash(:,startFrom(1):startFrom(1)+3,i),2);
            datarun(walk).white.nd6=[mean(white_flash(:,startFrom(2):startFrom(2)+3,i),2) mean(white_flash(:,startFrom(3):startFrom(3)+3,i),2) mean(white_flash(:,startFrom(4):startFrom(4)+3,i),2)];
            datarun(walk).white.nd5=[mean(white_flash(:,startFrom(5):startFrom(5)+3,i),2) mean(white_flash(:,startFrom(6):startFrom(6)+3,i),2) mean(white_flash(:,startFrom(7):startFrom(7)+3,i),2)];
            datarun(walk).white.nd4=[mean(white_flash(:,startFrom(8):startFrom(8)+3,i),2) mean(white_flash(:,startFrom(9):startFrom(9)+3,i),2) mean(white_flash(:,startFrom(10):startFrom(10)+3,i),2)];

            datarun(walk).chirp.nd7=mean(chirp(:,startFrom(1):startFrom(1)+3,i),2);
            datarun(walk).chirp.nd6=[mean(chirp(:,startFrom(2):startFrom(2)+3,i),2) mean(chirp(:,startFrom(3):startFrom(3)+3,i),2) mean(chirp(:,startFrom(4):startFrom(4)+3,i),2)];
            datarun(walk).chirp.nd5=[mean(chirp(:,startFrom(5):startFrom(5)+3,i),2) mean(chirp(:,startFrom(6):startFrom(6)+3,i),2) mean(chirp(:,startFrom(7):startFrom(7)+3,i),2)];
            datarun(walk).chirp.nd4=[mean(chirp(:,startFrom(8):startFrom(8)+3,i),2) mean(chirp(:,startFrom(9):startFrom(9)+3,i),2) mean(chirp(:,startFrom(10):startFrom(10)+3,i),2)];

            
            datarun(walk).lf.nd7=mean(LinearFilter(:,startFrom(1):startFrom(1)+3,i),2);
            datarun(walk).lf.nd6=[mean(LinearFilter(:,startFrom(2):startFrom(2)+3,i),2) mean(LinearFilter(:,startFrom(3):startFrom(3)+3,i),2) mean(LinearFilter(:,startFrom(4):startFrom(4)+3,i),2)];
            datarun(walk).lf.nd5=[mean(LinearFilter(:,startFrom(5):startFrom(5)+3,i),2) mean(LinearFilter(:,startFrom(6):startFrom(6)+3,i),2) mean(LinearFilter(:,startFrom(7):startFrom(7)+3,i),2)];
            datarun(walk).lf.nd4=[mean(LinearFilter(:,startFrom(8):startFrom(8)+3,i),2) mean(LinearFilter(:,startFrom(9):startFrom(9)+3,i),2) mean(LinearFilter(:,startFrom(10):startFrom(10)+3,i),2)];
           
            
            walk=walk+1;
        end
    end    
       
end

datarun_nd=datarun;
save('/Users/alexth/Desktop/old_stuff/ONresp/datarun_nd.mat','datarun_nd')


%% Cumulative data - averaged (checked for discarded units) by nd PLUS CONSISTENCY CHECK

clear

dates{1}='2013-09-03';
dates{2}='20130723a';
dates{3}='20130723b';
dates{4}='20130906';

datarun=struct('date',{},'name',{},'white',{},'black',{},'lf',{},'chirp',{},'polarity',{},'consistency',{});

walk=1;
for ddd=1:4
    
    date=dates{ddd};
    path2load=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
    load([path2load,'quick'])
    load([path2load,'LF_HC'])
    load([path2load,'chirp'])
    
    path2check=['/Users/alexth/Desktop/old_stuff/ONresp/cplot/'];
    names2check=dir(path2check);
    
    if strcmp(date,'20130723a')
        startFrom=[3 7 15 24 28 36 47 51 59 71];
    elseif strcmp(date,'2013-09-03')
        startFrom=[3 7 15 27 31 39 51 55 63 77];
    else
        startFrom=[3 7 15 27 31 39 51 55 63 75]; % number of trial to start averaging. First is for ND7, then 3 for ND6 (control, APB, wash), 3 fopr ND5 and 3 for ND4. Take 4 trials for each average from the middle of the run.
    end
    
    for i=1:size(black_flash,3)
        
        j=1;
        while ~strcmp(names{i},names2check(j).name(1:end-4)) && j<length(names2check)
            j=j+1;
        end        
        if strcmp(names{i},names2check(j).name(1:end-4))
            datarun(walk).date=date;
            datarun(walk).name=names{i};
            
            datarun(walk).black.nd7=mean(black_flash(:,startFrom(1):startFrom(1)+3,i),2);
            datarun(walk).black.nd6=[mean(black_flash(:,startFrom(2):startFrom(2)+3,i),2) mean(black_flash(:,startFrom(3):startFrom(3)+3,i),2) mean(black_flash(:,startFrom(4):startFrom(4)+3,i),2)];
            datarun(walk).black.nd5=[mean(black_flash(:,startFrom(5):startFrom(5)+3,i),2) mean(black_flash(:,startFrom(6):startFrom(6)+3,i),2) mean(black_flash(:,startFrom(7):startFrom(7)+3,i),2)];
            datarun(walk).black.nd4=[mean(black_flash(:,startFrom(8):startFrom(8)+3,i),2) mean(black_flash(:,startFrom(9):startFrom(9)+3,i),2) mean(black_flash(:,startFrom(10):startFrom(10)+3,i),2)];
            datarun(walk).white.nd7=mean(white_flash(:,startFrom(1):startFrom(1)+3,i),2);
            datarun(walk).white.nd6=[mean(white_flash(:,startFrom(2):startFrom(2)+3,i),2) mean(white_flash(:,startFrom(3):startFrom(3)+3,i),2) mean(white_flash(:,startFrom(4):startFrom(4)+3,i),2)];
            datarun(walk).white.nd5=[mean(white_flash(:,startFrom(5):startFrom(5)+3,i),2) mean(white_flash(:,startFrom(6):startFrom(6)+3,i),2) mean(white_flash(:,startFrom(7):startFrom(7)+3,i),2)];
            datarun(walk).white.nd4=[mean(white_flash(:,startFrom(8):startFrom(8)+3,i),2) mean(white_flash(:,startFrom(9):startFrom(9)+3,i),2) mean(white_flash(:,startFrom(10):startFrom(10)+3,i),2)];

            datarun(walk).chirp.nd7=mean(chirp(:,startFrom(1):startFrom(1)+3,i),2);
            datarun(walk).chirp.nd6=[mean(chirp(:,startFrom(2):startFrom(2)+3,i),2) mean(chirp(:,startFrom(3):startFrom(3)+3,i),2) mean(chirp(:,startFrom(4):startFrom(4)+3,i),2)];
            datarun(walk).chirp.nd5=[mean(chirp(:,startFrom(5):startFrom(5)+3,i),2) mean(chirp(:,startFrom(6):startFrom(6)+3,i),2) mean(chirp(:,startFrom(7):startFrom(7)+3,i),2)];
            datarun(walk).chirp.nd4=[mean(chirp(:,startFrom(8):startFrom(8)+3,i),2) mean(chirp(:,startFrom(9):startFrom(9)+3,i),2) mean(chirp(:,startFrom(10):startFrom(10)+3,i),2)];

            
            datarun(walk).lf.nd7=mean(LinearFilter(:,startFrom(1):startFrom(1)+3,i),2);
            datarun(walk).lf.nd6=[mean(LinearFilter(:,startFrom(2):startFrom(2)+3,i),2) mean(LinearFilter(:,startFrom(3):startFrom(3)+3,i),2) mean(LinearFilter(:,startFrom(4):startFrom(4)+3,i),2)];
            datarun(walk).lf.nd5=[mean(LinearFilter(:,startFrom(5):startFrom(5)+3,i),2) mean(LinearFilter(:,startFrom(6):startFrom(6)+3,i),2) mean(LinearFilter(:,startFrom(7):startFrom(7)+3,i),2)];
            datarun(walk).lf.nd4=[mean(LinearFilter(:,startFrom(8):startFrom(8)+3,i),2) mean(LinearFilter(:,startFrom(9):startFrom(9)+3,i),2) mean(LinearFilter(:,startFrom(10):startFrom(10)+3,i),2)];
  

            % consistency measure
            datarun(walk).consistency.black.nd7=consist_check(black_flash(:,startFrom(1):startFrom(1)+3,i));
            datarun(walk).consistency.black.nd6=[consist_check(black_flash(:,startFrom(2):startFrom(2)+3,i)) consist_check(black_flash(:,startFrom(3):startFrom(3)+3,i)) consist_check(black_flash(:,startFrom(4):startFrom(4)+3,i))];
            datarun(walk).consistency.black.nd5=[consist_check(black_flash(:,startFrom(5):startFrom(5)+3,i)) consist_check(black_flash(:,startFrom(6):startFrom(6)+3,i)) consist_check(black_flash(:,startFrom(7):startFrom(7)+3,i))];
            datarun(walk).consistency.black.nd4=[consist_check(black_flash(:,startFrom(8):startFrom(8)+3,i)) consist_check(black_flash(:,startFrom(9):startFrom(9)+3,i)) consist_check(black_flash(:,startFrom(10):startFrom(10)+3,i))];
            datarun(walk).consistency.white.nd7=consist_check(white_flash(:,startFrom(1):startFrom(1)+3,i));
            datarun(walk).consistency.white.nd6=[consist_check(white_flash(:,startFrom(2):startFrom(2)+3,i)) consist_check(white_flash(:,startFrom(3):startFrom(3)+3,i)) consist_check(white_flash(:,startFrom(4):startFrom(4)+3,i))];
            datarun(walk).consistency.white.nd5=[consist_check(white_flash(:,startFrom(5):startFrom(5)+3,i)) consist_check(white_flash(:,startFrom(6):startFrom(6)+3,i)) consist_check(white_flash(:,startFrom(7):startFrom(7)+3,i))];
            datarun(walk).consistency.white.nd4=[consist_check(white_flash(:,startFrom(8):startFrom(8)+3,i)) consist_check(white_flash(:,startFrom(9):startFrom(9)+3,i)) consist_check(white_flash(:,startFrom(10):startFrom(10)+3,i))];

            datarun(walk).consistency.chirp.nd7=consist_check(chirp(:,startFrom(1):startFrom(1)+3,i));
            datarun(walk).consistency.chirp.nd6=[consist_check(chirp(:,startFrom(2):startFrom(2)+3,i)) consist_check(chirp(:,startFrom(3):startFrom(3)+3,i)) consist_check(chirp(:,startFrom(4):startFrom(4)+3,i))];
            datarun(walk).consistency.chirp.nd5=[consist_check(chirp(:,startFrom(5):startFrom(5)+3,i)) consist_check(chirp(:,startFrom(6):startFrom(6)+3,i)) consist_check(chirp(:,startFrom(7):startFrom(7)+3,i))];
            datarun(walk).consistency.chirp.nd4=[consist_check(chirp(:,startFrom(8):startFrom(8)+3,i)) consist_check(chirp(:,startFrom(9):startFrom(9)+3,i)) consist_check(chirp(:,startFrom(10):startFrom(10)+3,i))];

            
            datarun(walk).consistency.lf.nd7=consist_check(LinearFilter(:,startFrom(1):startFrom(1)+3,i));
            datarun(walk).consistency.lf.nd6=[consist_check(LinearFilter(:,startFrom(2):startFrom(2)+3,i)) consist_check(LinearFilter(:,startFrom(3):startFrom(3)+3,i)) consist_check(LinearFilter(:,startFrom(4):startFrom(4)+3,i))];
            datarun(walk).consistency.lf.nd5=[consist_check(LinearFilter(:,startFrom(5):startFrom(5)+3,i)) consist_check(LinearFilter(:,startFrom(6):startFrom(6)+3,i)) consist_check(LinearFilter(:,startFrom(7):startFrom(7)+3,i))];
            datarun(walk).consistency.lf.nd4=[consist_check(LinearFilter(:,startFrom(8):startFrom(8)+3,i)) consist_check(LinearFilter(:,startFrom(9):startFrom(9)+3,i)) consist_check(LinearFilter(:,startFrom(10):startFrom(10)+3,i))];

            
            walk=walk+1;
        end
    end    
       
end

tmp=zeros(500,4,size(datarun,2));

for i=1:size(datarun,2)
    k=[datarun(i).lf.nd7 datarun(i).lf.nd6(:,1) datarun(i).lf.nd5(:,1) datarun(i).lf.nd4(:,1)];
    k=k-repmat(mean(k),500,1);
    k=k./abs(repmat(max(k)-min(k),500,1));
    tmp(:,:,i)=k;
end

tmp(isnan(tmp))=0;

for i=1:4
    t=squeeze(tmp(50:250,i,:));
    x=t';
    x=cov(x);
    [V,~]=eig(x);
    pc_vectors=V(:,end);
    pc1(:,i)=t'*pc_vectors;
end


onCells=find(sum(pc1<0,2)>2); % ON cells
offCells=find(sum(pc1>0,2)>2); % OFF cells

strangeCells=find(sum(pc1>0,2)==2); % changing cells
clear tmp
tmp(onCells)=1;
tmp(offCells)=-1;
tmp(strangeCells)=0;

for i=1:size(datarun,2)
    datarun(i).polarity=tmp(i);
end

datarun_nd=datarun;

save('/Users/alexth/Desktop/old_stuff/ONresp/datarun_nd.mat','datarun_nd')


%% Black flash OFF cells
clear
load('/Users/alexth/Desktop/old_stuff/ONresp/datarun_nd.mat')

for i=1:49
    subplot(7,7,i)
    plot(datarun_nd(i).black.nd6)
end

clear offBlack onDelBlack

for i=1:111
    tmp=datarun_nd(i).black.nd6;
    offBlack(i,1:3)=max(tmp(501:1000,:))-max(tmp(1:500,:));
end
a=offBlack(:,1)>5;
sum(a)
b=find(a);
totalCells=length(b);

for i=1:sum(a)
    subplot(8,8,i)
    plot(datarun_nd(b(i)).black.nd6)
end
    
    
for i=1:111
    tmp=datarun_nd(i).black.nd6;
    offBlack(i,1:3)=max(tmp(501:1000,:))-max(tmp(1:500,:));
end

for i=1:111
    tmp=datarun_nd(i).black.nd6;
    onDelBlack(i,1:3)=max(tmp(3250:3600,:))-max(tmp(1:500,:));
end

for i=1:111
    tmp=datarun_nd(i).black.nd6;
    onEarBlack(i,1:3)=max(tmp(2600:2800,:))-max(tmp(1:500,:));
end

figure

templateResp=onDelBlack;
startSbpl=1;

tmp1=templateResp(b,:);
[~,c]=sort(tmp1(:,1));

subplot(3,3,startSbpl)
plot(templateResp(b(c),:))
legend('Control','APB','wash')
line([0 totalCells],[0 0],'color','k')
line([0 totalCells],[5 5],'LineStyle','--','color','k')
text(50,7,'threshold')
title([{'Delayed ON response'},{'sorted by control amplitude'}])
xlabel('cell ordinal number')
ylabel('amplitude of response, Hz')
axis([0 totalCells -30 100])

subplot(3,3,startSbpl+1)
bar(sum(templateResp(b,:)>5))
set(gca,'xticklabel',{'Control','APB','wash'})
xlabel('N of cells')
title('Delayed ON response >5Hz')
axis([0 4 0 totalCells])

subplot(3,3,startSbpl+2)
amplCh(:,1)=templateResp(b,1)-templateResp(b,2);
amplCh(:,2)=templateResp(b,1)-templateResp(b,3);
plot(templateResp(b,1),amplCh(:,1),'*');
hold on
plot(templateResp(b,1),amplCh(:,2),'*g');
line([-100 100],[0 0],'color','k')
line([0 0],[-100 100],'color','k')
legend('Control-APB','Control-wash','location','best')
xlabel('Control Amplitude')
ylabel('Amplitude difference')
title('Delayed ON response amplitude change')



templateResp=onEarBlack;
startSbpl=4;

tmp1=templateResp(b,:);
[~,c]=sort(tmp1(:,1));

subplot(3,3,startSbpl)
plot(templateResp(b(c),:))
legend('Control','APB','wash')
line([0 totalCells],[0 0],'color','k')
line([0 totalCells],[5 5],'LineStyle','--','color','k')
text(50,7,'threshold')
title([{'Early ON response'},{'sorted by control amplitude'}])
xlabel('cell ordinal number')
ylabel('amplitude of response, Hz')
axis([0 totalCells -30 100])

subplot(3,3,startSbpl+1)
bar(sum(templateResp(b,:)>5))
set(gca,'xticklabel',{'Control','APB','wash'})
xlabel('N of cells')
title('Early ON response >5Hz')
axis([0 4 0 totalCells])

subplot(3,3,startSbpl+2)
amplCh(:,1)=templateResp(b,1)-templateResp(b,2);
amplCh(:,2)=templateResp(b,1)-templateResp(b,3);
plot(templateResp(b,1),amplCh(:,1),'*');
hold on
plot(templateResp(b,1),amplCh(:,2),'*g');
line([-100 100],[0 0],'color','k')
line([0 0],[-100 100],'color','k')
legend('Control-APB','Control-wash','location','best')
xlabel('Control Amplitude')
ylabel('Amplitude difference')
title('Early ON response amplitude change')




templateResp=offBlack;
startSbpl=7;

tmp1=templateResp(b,:);
[~,c]=sort(tmp1(:,1));

subplot(3,3,startSbpl)
plot(templateResp(b(c),:))
legend('Control','APB','wash')
line([0 totalCells],[0 0],'color','k')
line([0 totalCells],[5 5],'LineStyle','--','color','k')
text(50,7,'threshold')
title([{'OFF response'},{'sorted by control amplitude'}])
xlabel('cell ordinal number')
ylabel('amplitude of response, Hz')
axis([0 totalCells -30 100])

subplot(3,3,startSbpl+1)
bar(sum(templateResp(b,:)>5))
set(gca,'xticklabel',{'Control','APB','wash'})
xlabel('N of cells')
title('OFF response >5Hz')
axis([0 4 0 totalCells])

subplot(3,3,startSbpl+2)
amplCh(:,1)=templateResp(b,1)-templateResp(b,2);
amplCh(:,2)=templateResp(b,1)-templateResp(b,3);
plot(templateResp(b,1),amplCh(:,1),'*');
hold on
plot(templateResp(b,1),amplCh(:,2),'*g');
line([-100 100],[0 0],'color','k')
line([0 0],[-100 100],'color','k')
legend('Control-APB','Control-wash','location','best')
xlabel('Control Amplitude')
ylabel('Amplitude difference')
title('OFF response amplitude change')






%% White flash OFF cells
clear
load('/Users/alexth/Desktop/old_stuff/ONresp/datarun_nd.mat')

for i=1:49
    subplot(7,7,i)
    plot(datarun_nd(i).white.nd6)
end

clear offBlack onDelBlack

for i=1:111
    tmp=datarun_nd(i).white.nd6;
    offWhite(i,1:3)=max(tmp(2501:3000,:))-max(tmp(1:500,:));
end
a=offWhite(:,1)>5;
sum(a)
b=find(a);
totalCells=length(b);

for i=1:sum(a)
    subplot(8,8,i)
    plot(datarun_nd(b(i)).white.nd6)
end
    
    
for i=1:111
    tmp=datarun_nd(i).white.nd6;
    offWhite(i,1:3)=max(tmp(2501:3000,:))-max(tmp(1:500,:));
end

for i=1:111
    tmp=datarun_nd(i).white.nd6;
    onDelWhite(i,1:3)=max(tmp(1250:1600,:))-max(tmp(1:500,:));
end

for i=1:111
    tmp=datarun_nd(i).white.nd6;
    onEarWhite(i,1:3)=max(tmp(600:800,:))-max(tmp(1:500,:));
end

figure

templateResp=onDelWhite;
startSbpl=1;

tmp1=templateResp(b,:);
[~,c]=sort(tmp1(:,1));

subplot(3,3,startSbpl)
plot(templateResp(b(c),:))
legend('Control','APB','wash')
line([0 totalCells],[0 0],'color','k')
line([0 totalCells],[5 5],'LineStyle','--','color','k')
text(50,7,'threshold')
title([{'Delayed ON response'},{'sorted by control amplitude'}])
xlabel('cell ordinal number')
ylabel('amplitude of response, Hz')
axis([0 totalCells -30 100])

subplot(3,3,startSbpl+1)
bar(sum(templateResp(b,:)>5))
set(gca,'xticklabel',{'Control','APB','wash'})
xlabel('N of cells')
title('Delayed ON response >5Hz')
axis([0 4 0 totalCells])

subplot(3,3,startSbpl+2)
amplCh(:,1)=templateResp(b,1)-templateResp(b,2);
amplCh(:,2)=templateResp(b,1)-templateResp(b,3);
plot(templateResp(b,1),amplCh(:,1),'*');
hold on
plot(templateResp(b,1),amplCh(:,2),'*g');
line([-100 100],[0 0],'color','k')
line([0 0],[-100 100],'color','k')
legend('Control-APB','Control-wash','location','best')
xlabel('Control Amplitude')
ylabel('Amplitude difference')
title('Delayed ON response amplitude change')



templateResp=onEarWhite;
startSbpl=4;

tmp1=templateResp(b,:);
[~,c]=sort(tmp1(:,1));

subplot(3,3,startSbpl)
plot(templateResp(b(c),:))
legend('Control','APB','wash')
line([0 totalCells],[0 0],'color','k')
line([0 totalCells],[5 5],'LineStyle','--','color','k')
text(50,7,'threshold')
title([{'Early ON response'},{'sorted by control amplitude'}])
xlabel('cell ordinal number')
ylabel('amplitude of response, Hz')
axis([0 totalCells -30 100])

subplot(3,3,startSbpl+1)
bar(sum(templateResp(b,:)>5))
set(gca,'xticklabel',{'Control','APB','wash'})
xlabel('N of cells')
title('Early ON response >5Hz')
axis([0 4 0 totalCells])

subplot(3,3,startSbpl+2)
amplCh(:,1)=templateResp(b,1)-templateResp(b,2);
amplCh(:,2)=templateResp(b,1)-templateResp(b,3);
plot(templateResp(b,1),amplCh(:,1),'*');
hold on
plot(templateResp(b,1),amplCh(:,2),'*g');
line([-100 100],[0 0],'color','k')
line([0 0],[-100 100],'color','k')
legend('Control-APB','Control-wash','location','best')
xlabel('Control Amplitude')
ylabel('Amplitude difference')
title('Early ON response amplitude change')




templateResp=offWhite;
startSbpl=7;

tmp1=templateResp(b,:);
[~,c]=sort(tmp1(:,1));

subplot(3,3,startSbpl)
plot(templateResp(b(c),:))
legend('Control','APB','wash')
line([0 totalCells],[0 0],'color','k')
line([0 totalCells],[5 5],'LineStyle','--','color','k')
text(50,7,'threshold')
title([{'OFF response'},{'sorted by control amplitude'}])
xlabel('cell ordinal number')
ylabel('amplitude of response, Hz')
axis([0 totalCells -30 100])

subplot(3,3,startSbpl+1)
bar(sum(templateResp(b,:)>5))
set(gca,'xticklabel',{'Control','APB','wash'})
xlabel('N of cells')
title('OFF response >5Hz')
axis([0 4 0 totalCells])

subplot(3,3,startSbpl+2)
amplCh(:,1)=templateResp(b,1)-templateResp(b,2);
amplCh(:,2)=templateResp(b,1)-templateResp(b,3);
plot(templateResp(b,1),amplCh(:,1),'*');
hold on
plot(templateResp(b,1),amplCh(:,2),'*g');
line([-100 100],[0 0],'color','k')
line([0 0],[-100 100],'color','k')
legend('Control-APB','Control-wash','location','best')
xlabel('Control Amplitude')
ylabel('Amplitude difference')
title('OFF response amplitude change')





%% White flash ON cells
clear
load('/Users/alexth/Desktop/old_stuff/ONresp/datarun_nd.mat')

for i=1:49
    subplot(7,7,i)
    plot(datarun_nd(i).white.nd6)
end


for i=1:111
    tmp=datarun_nd(i).white.nd6;
    onWhite(i,1:3)=max(tmp(501:1000,:))-max(tmp(1:500,:));
end
a=onWhite(:,1)>5;
sum(a)
b=find(a);
totalCells=length(b);

for i=1:sum(a)
    subplot(8,8,i)
    plot(datarun_nd(b(i)).white.nd6)
end
    
    
for i=1:111
    tmp=datarun_nd(i).white.nd6;
    onWhite(i,1:3)=max(tmp(501:1000,:))-max(tmp(1:500,:));
end

figure

rowRespWhite=0;
rowRespBlack=0;
for i=1:totalCells
    tmp=datarun_nd(b(i)).white.nd6;
    rowRespWhite=rowRespWhite+tmp;
    tmp=datarun_nd(b(i)).black.nd6;
    rowRespBlack=rowRespBlack+tmp;
end

subplot(1,2,1)
plot(rowRespWhite/totalCells)
subplot(1,2,2)
plot(rowRespBlack/totalCells)




%% OFF cells early ON response compare at ND7
clear
load('/Users/alexth/Desktop/old_stuff/ONresp/datarun_nd.mat')

for i=1:49
    subplot(7,7,i)
    plot(datarun_nd(i).black.nd6)
end

clear offBlack onDelBlack

for i=1:111
    tmp=datarun_nd(i).black.nd6;
    offBlack(i,1:3)=max(tmp(501:1000,:))-max(tmp(1:500,:));
end
a=offBlack(:,1)>5;
sum(a)
b=find(a);
totalCells=length(b);

for i=1:sum(a)
    subplot(8,8,i)
    plot(datarun_nd(b(i)).black.nd6)
end
    

for i=1:111
    tmp=datarun_nd(i).black.nd6;
    onEarBlack(i,1:3)=max(tmp(2600:2800,:))-max(tmp(1:500,:));
end

for i=1:111
    tmp=datarun_nd(i).black.nd7;
    onEarBlackND7(i)=max(tmp(2600:2800,:))-max(tmp(1:500,:));
end

for i=1:111
    tmp=datarun_nd(i).black.nd5;
    onEarBlackND5(i,1:3)=max(tmp(2550:2750,:))-max(tmp(1:500,:));
end

for i=1:111
    tmp=datarun_nd(i).black.nd4;
    onEarBlackND4(i,1:3)=max(tmp(2550:2750,:))-max(tmp(1:500,:));
end

% find 'true' ON-OFF cells - OFF with early ON responses at all NDs (control conditions)
trueOnOff=intersect(find(onEarBlack(:,1)>5&onEarBlackND7(:,1)>5&onEarBlackND5(:,1)>5&onEarBlackND4(:,1)>5),b);

% find 'pure' OFF cells with no ON responses at all NDs (control conditions)
pureOff=intersect(find(onEarBlack(:,1)<5&onEarBlackND7(:,1)<5&onEarBlackND5(:,1)<5&onEarBlackND4(:,1)<5),b);

% find 'putative' ON-OFF cells with ON responses at some NDs (control conditions)
putOnOff=setdiff(setdiff(b,trueOnOff),pureOff);


%    ######### Plot Figure ##########

%plot ND7
rowRespTrueOnOff=0;
for i=1:length(trueOnOff)
    rowRespTrueOnOff=rowRespTrueOnOff+[datarun_nd(trueOnOff(i)).black.nd7;...
       zeros(200,1); datarun_nd(trueOnOff(i)).white.nd7];
end
rowRespTrueOnOff=rowRespTrueOnOff/length(trueOnOff);

rowRespPureOff=0;
for i=1:length(pureOff)
    rowRespPureOff=rowRespPureOff+[datarun_nd(pureOff(i)).black.nd7;...
       zeros(200,1); datarun_nd(pureOff(i)).white.nd7];
end
rowRespPureOff=rowRespPureOff/length(pureOff);

rowRespPutOnOff=0;
for i=1:length(putOnOff)
    rowRespPutOnOff=rowRespPutOnOff+[datarun_nd(putOnOff(i)).black.nd7;...
       zeros(200,1); datarun_nd(putOnOff(i)).white.nd7];
end
rowRespPutOnOff=rowRespPutOnOff/length(putOnOff);

figure
subplot(4,3,1)
plot(rowRespTrueOnOff)
a=rowRespTrueOnOff;
%plot black flash timing
line([0 500], [max(a)*1.2, max(a)*1.2],'color','k')
line([500 500], [max(a)*1.1, max(a)*1.2],'color','k')
line([500 2500], [max(a)*1.1, max(a)*1.1],'color','k')
line([2500 2500], [max(a)*1.1, max(a)*1.2],'color','k')
line([2500 4500], [max(a)*1.2, max(a)*1.2],'color','k')
%plot white flash timing
line([4700 5300], [max(a)*1.2, max(a)*1.2],'color','k')
line([5300 5300], [max(a)*1.2, max(a)*1.3],'color','k')
line([5300 7300], [max(a)*1.3, max(a)*1.3],'color','k')
line([7300 7300], [max(a)*1.2, max(a)*1.3],'color','k')
line([7300 9300], [max(a)*1.2, max(a)*1.2],'color','k')
%fill gap
for i=1:15
    line([4500 4700], [max(a)*1.4/15*(i-1),max(a)*1.4/15*(i)],'color','k')
end
axis([0 9200 0 max(max(a)*1.35,1)])
title('ND7 true ON-OFF')

subplot(4,3,2)
plot(rowRespPureOff)
a=rowRespPureOff;
%plot black flash timing
line([0 500], [max(a)*1.2, max(a)*1.2],'color','k')
line([500 500], [max(a)*1.1, max(a)*1.2],'color','k')
line([500 2500], [max(a)*1.1, max(a)*1.1],'color','k')
line([2500 2500], [max(a)*1.1, max(a)*1.2],'color','k')
line([2500 4500], [max(a)*1.2, max(a)*1.2],'color','k')
%plot white flash timing
line([4700 5300], [max(a)*1.2, max(a)*1.2],'color','k')
line([5300 5300], [max(a)*1.2, max(a)*1.3],'color','k')
line([5300 7300], [max(a)*1.3, max(a)*1.3],'color','k')
line([7300 7300], [max(a)*1.2, max(a)*1.3],'color','k')
line([7300 9300], [max(a)*1.2, max(a)*1.2],'color','k')
%fill gap
for i=1:15
    line([4500 4700], [max(a)*1.4/15*(i-1),max(a)*1.4/15*(i)],'color','k')
end
axis([0 9200 0 max(max(a)*1.35,1)])
title('ND7 pure OFF')


subplot(4,3,3)
plot(rowRespPutOnOff)
a=rowRespPutOnOff;
%plot black flash timing
line([0 500], [max(a)*1.2, max(a)*1.2],'color','k')
line([500 500], [max(a)*1.1, max(a)*1.2],'color','k')
line([500 2500], [max(a)*1.1, max(a)*1.1],'color','k')
line([2500 2500], [max(a)*1.1, max(a)*1.2],'color','k')
line([2500 4500], [max(a)*1.2, max(a)*1.2],'color','k')
%plot white flash timing
line([4700 5300], [max(a)*1.2, max(a)*1.2],'color','k')
line([5300 5300], [max(a)*1.2, max(a)*1.3],'color','k')
line([5300 7300], [max(a)*1.3, max(a)*1.3],'color','k')
line([7300 7300], [max(a)*1.2, max(a)*1.3],'color','k')
line([7300 9300], [max(a)*1.2, max(a)*1.2],'color','k')
%fill gap
for i=1:15
    line([4500 4700], [max(a)*1.4/15*(i-1),max(a)*1.4/15*(i)],'color','k')
end
axis([0 9200 0 max(max(a)*1.35,1)])
title('ND7 putative OFF')
legend('control')


%plot ND6
rowRespTrueOnOff=0;
for i=1:length(trueOnOff)
    rowRespTrueOnOff=rowRespTrueOnOff+[datarun_nd(trueOnOff(i)).black.nd6;...
       zeros(200,3); datarun_nd(trueOnOff(i)).white.nd6];
end
rowRespTrueOnOff=rowRespTrueOnOff/length(trueOnOff);

rowRespPureOff=0;
for i=1:length(pureOff)
    rowRespPureOff=rowRespPureOff+[datarun_nd(pureOff(i)).black.nd6;...
       zeros(200,3); datarun_nd(pureOff(i)).white.nd6];
end
rowRespPureOff=rowRespPureOff/length(pureOff);

rowRespPutOnOff=0;
for i=1:length(putOnOff)
    rowRespPutOnOff=rowRespPutOnOff+[datarun_nd(putOnOff(i)).black.nd6;...
       zeros(200,3); datarun_nd(putOnOff(i)).white.nd6];
end
rowRespPutOnOff=rowRespPutOnOff/length(putOnOff);


subplot(4,3,4)
plot(rowRespTrueOnOff)
a=max(rowRespTrueOnOff');
%plot black flash timing
line([0 500], [max(a)*1.2, max(a)*1.2],'color','k')
line([500 500], [max(a)*1.1, max(a)*1.2],'color','k')
line([500 2500], [max(a)*1.1, max(a)*1.1],'color','k')
line([2500 2500], [max(a)*1.1, max(a)*1.2],'color','k')
line([2500 4500], [max(a)*1.2, max(a)*1.2],'color','k')
%plot white flash timing
line([4700 5300], [max(a)*1.2, max(a)*1.2],'color','k')
line([5300 5300], [max(a)*1.2, max(a)*1.3],'color','k')
line([5300 7300], [max(a)*1.3, max(a)*1.3],'color','k')
line([7300 7300], [max(a)*1.2, max(a)*1.3],'color','k')
line([7300 9300], [max(a)*1.2, max(a)*1.2],'color','k')
%fill gap
for i=1:15
    line([4500 4700], [max(a)*1.4/15*(i-1),max(a)*1.4/15*(i)],'color','k')
end
axis([0 9200 0 max(max(a)*1.35,1)])
title('ND6 true ON-OFF')

subplot(4,3,5)
plot(rowRespPureOff)
a=max(rowRespPureOff');
%plot black flash timing
line([0 500], [max(a)*1.2, max(a)*1.2],'color','k')
line([500 500], [max(a)*1.1, max(a)*1.2],'color','k')
line([500 2500], [max(a)*1.1, max(a)*1.1],'color','k')
line([2500 2500], [max(a)*1.1, max(a)*1.2],'color','k')
line([2500 4500], [max(a)*1.2, max(a)*1.2],'color','k')
%plot white flash timing
line([4700 5300], [max(a)*1.2, max(a)*1.2],'color','k')
line([5300 5300], [max(a)*1.2, max(a)*1.3],'color','k')
line([5300 7300], [max(a)*1.3, max(a)*1.3],'color','k')
line([7300 7300], [max(a)*1.2, max(a)*1.3],'color','k')
line([7300 9300], [max(a)*1.2, max(a)*1.2],'color','k')
%fill gap
for i=1:15
    line([4500 4700], [max(a)*1.4/15*(i-1),max(a)*1.4/15*(i)],'color','k')
end
axis([0 9200 0 max(max(a)*1.35,1)])
title('ND6 pure OFF')


subplot(4,3,6)
plot(rowRespPutOnOff)
legend('control','apb','wash')
a=max(rowRespPutOnOff');
%plot black flash timing
line([0 500], [max(a)*1.2, max(a)*1.2],'color','k')
line([500 500], [max(a)*1.1, max(a)*1.2],'color','k')
line([500 2500], [max(a)*1.1, max(a)*1.1],'color','k')
line([2500 2500], [max(a)*1.1, max(a)*1.2],'color','k')
line([2500 4500], [max(a)*1.2, max(a)*1.2],'color','k')
%plot white flash timing
line([4700 5300], [max(a)*1.2, max(a)*1.2],'color','k')
line([5300 5300], [max(a)*1.2, max(a)*1.3],'color','k')
line([5300 7300], [max(a)*1.3, max(a)*1.3],'color','k')
line([7300 7300], [max(a)*1.2, max(a)*1.3],'color','k')
line([7300 9300], [max(a)*1.2, max(a)*1.2],'color','k')
%fill gap
for i=1:15
    line([4500 4700], [max(a)*1.4/15*(i-1),max(a)*1.4/15*(i)],'color','k')
end
axis([0 9200 0 max(max(a)*1.35,1)])
title('ND6 putative OFF')




%plot ND5
rowRespTrueOnOff=0;
for i=1:length(trueOnOff)
    rowRespTrueOnOff=rowRespTrueOnOff+[datarun_nd(trueOnOff(i)).black.nd5;...
       zeros(200,3); datarun_nd(trueOnOff(i)).white.nd5];
end
rowRespTrueOnOff=rowRespTrueOnOff/length(trueOnOff);

rowRespPureOff=0;
for i=1:length(pureOff)
    rowRespPureOff=rowRespPureOff+[datarun_nd(pureOff(i)).black.nd5;...
       zeros(200,3); datarun_nd(pureOff(i)).white.nd5];
end
rowRespPureOff=rowRespPureOff/length(pureOff);

rowRespPutOnOff=0;
for i=1:length(putOnOff)
    rowRespPutOnOff=rowRespPutOnOff+[datarun_nd(putOnOff(i)).black.nd5;...
       zeros(200,3); datarun_nd(putOnOff(i)).white.nd5];
end
rowRespPutOnOff=rowRespPutOnOff/length(putOnOff);


subplot(4,3,7)
plot(rowRespTrueOnOff)
a=max(rowRespTrueOnOff');
%plot black flash timing
line([0 500], [max(a)*1.2, max(a)*1.2],'color','k')
line([500 500], [max(a)*1.1, max(a)*1.2],'color','k')
line([500 2500], [max(a)*1.1, max(a)*1.1],'color','k')
line([2500 2500], [max(a)*1.1, max(a)*1.2],'color','k')
line([2500 4500], [max(a)*1.2, max(a)*1.2],'color','k')
%plot white flash timing
line([4700 5300], [max(a)*1.2, max(a)*1.2],'color','k')
line([5300 5300], [max(a)*1.2, max(a)*1.3],'color','k')
line([5300 7300], [max(a)*1.3, max(a)*1.3],'color','k')
line([7300 7300], [max(a)*1.2, max(a)*1.3],'color','k')
line([7300 9300], [max(a)*1.2, max(a)*1.2],'color','k')
%fill gap
for i=1:15
    line([4500 4700], [max(a)*1.4/15*(i-1),max(a)*1.4/15*(i)],'color','k')
end
axis([0 9200 0 max(max(a)*1.35,1)])
title('ND5 true ON-OFF')

subplot(4,3,8)
plot(rowRespPureOff)
a=max(rowRespPureOff');
%plot black flash timing
line([0 500], [max(a)*1.2, max(a)*1.2],'color','k')
line([500 500], [max(a)*1.1, max(a)*1.2],'color','k')
line([500 2500], [max(a)*1.1, max(a)*1.1],'color','k')
line([2500 2500], [max(a)*1.1, max(a)*1.2],'color','k')
line([2500 4500], [max(a)*1.2, max(a)*1.2],'color','k')
%plot white flash timing
line([4700 5300], [max(a)*1.2, max(a)*1.2],'color','k')
line([5300 5300], [max(a)*1.2, max(a)*1.3],'color','k')
line([5300 7300], [max(a)*1.3, max(a)*1.3],'color','k')
line([7300 7300], [max(a)*1.2, max(a)*1.3],'color','k')
line([7300 9300], [max(a)*1.2, max(a)*1.2],'color','k')
%fill gap
for i=1:15
    line([4500 4700], [max(a)*1.4/15*(i-1),max(a)*1.4/15*(i)],'color','k')
end
axis([0 9200 0 max(max(a)*1.35,1)])
title('ND5 pure OFF')


subplot(4,3,9)
plot(rowRespPutOnOff)
legend('control','apb','wash')
a=max(rowRespPutOnOff');
%plot black flash timing
line([0 500], [max(a)*1.2, max(a)*1.2],'color','k')
line([500 500], [max(a)*1.1, max(a)*1.2],'color','k')
line([500 2500], [max(a)*1.1, max(a)*1.1],'color','k')
line([2500 2500], [max(a)*1.1, max(a)*1.2],'color','k')
line([2500 4500], [max(a)*1.2, max(a)*1.2],'color','k')
%plot white flash timing
line([4700 5300], [max(a)*1.2, max(a)*1.2],'color','k')
line([5300 5300], [max(a)*1.2, max(a)*1.3],'color','k')
line([5300 7300], [max(a)*1.3, max(a)*1.3],'color','k')
line([7300 7300], [max(a)*1.2, max(a)*1.3],'color','k')
line([7300 9300], [max(a)*1.2, max(a)*1.2],'color','k')
%fill gap
for i=1:15
    line([4500 4700], [max(a)*1.4/15*(i-1),max(a)*1.4/15*(i)],'color','k')
end
axis([0 9200 0 max(max(a)*1.35,1)])
title('ND5 putative OFF')




%plot ND4
rowRespTrueOnOff=0;
for i=1:length(trueOnOff)
    rowRespTrueOnOff=rowRespTrueOnOff+[datarun_nd(trueOnOff(i)).black.nd4;...
       zeros(200,3); datarun_nd(trueOnOff(i)).white.nd4];
end
rowRespTrueOnOff=rowRespTrueOnOff/length(trueOnOff);

rowRespPureOff=0;
for i=1:length(pureOff)
    rowRespPureOff=rowRespPureOff+[datarun_nd(pureOff(i)).black.nd4;...
       zeros(200,3); datarun_nd(pureOff(i)).white.nd4];
end
rowRespPureOff=rowRespPureOff/length(pureOff);

rowRespPutOnOff=0;
for i=1:length(putOnOff)
    rowRespPutOnOff=rowRespPutOnOff+[datarun_nd(putOnOff(i)).black.nd4;...
       zeros(200,3); datarun_nd(putOnOff(i)).white.nd4];
end
rowRespPutOnOff=rowRespPutOnOff/length(putOnOff);


subplot(4,3,10)
plot(rowRespTrueOnOff)
a=max(rowRespTrueOnOff');
%plot black flash timing
line([0 500], [max(a)*1.2, max(a)*1.2],'color','k')
line([500 500], [max(a)*1.1, max(a)*1.2],'color','k')
line([500 2500], [max(a)*1.1, max(a)*1.1],'color','k')
line([2500 2500], [max(a)*1.1, max(a)*1.2],'color','k')
line([2500 4500], [max(a)*1.2, max(a)*1.2],'color','k')
%plot white flash timing
line([4700 5300], [max(a)*1.2, max(a)*1.2],'color','k')
line([5300 5300], [max(a)*1.2, max(a)*1.3],'color','k')
line([5300 7300], [max(a)*1.3, max(a)*1.3],'color','k')
line([7300 7300], [max(a)*1.2, max(a)*1.3],'color','k')
line([7300 9300], [max(a)*1.2, max(a)*1.2],'color','k')
%fill gap
for i=1:15
    line([4500 4700], [max(a)*1.4/15*(i-1),max(a)*1.4/15*(i)],'color','k')
end
axis([0 9200 0 max(max(a)*1.35,1)])
title('ND4 true ON-OFF')

subplot(4,3,11)
plot(rowRespPureOff)
a=max(rowRespPureOff');
%plot black flash timing
line([0 500], [max(a)*1.2, max(a)*1.2],'color','k')
line([500 500], [max(a)*1.1, max(a)*1.2],'color','k')
line([500 2500], [max(a)*1.1, max(a)*1.1],'color','k')
line([2500 2500], [max(a)*1.1, max(a)*1.2],'color','k')
line([2500 4500], [max(a)*1.2, max(a)*1.2],'color','k')
%plot white flash timing
line([4700 5300], [max(a)*1.2, max(a)*1.2],'color','k')
line([5300 5300], [max(a)*1.2, max(a)*1.3],'color','k')
line([5300 7300], [max(a)*1.3, max(a)*1.3],'color','k')
line([7300 7300], [max(a)*1.2, max(a)*1.3],'color','k')
line([7300 9300], [max(a)*1.2, max(a)*1.2],'color','k')
%fill gap
for i=1:15
    line([4500 4700], [max(a)*1.4/15*(i-1),max(a)*1.4/15*(i)],'color','k')
end
axis([0 9200 0 max(max(a)*1.35,1)])
title('ND4 pure OFF')


subplot(4,3,12)
plot(rowRespPutOnOff)
legend('control','apb','wash')
a=max(rowRespPutOnOff');
%plot black flash timing
line([0 500], [max(a)*1.2, max(a)*1.2],'color','k')
line([500 500], [max(a)*1.1, max(a)*1.2],'color','k')
line([500 2500], [max(a)*1.1, max(a)*1.1],'color','k')
line([2500 2500], [max(a)*1.1, max(a)*1.2],'color','k')
line([2500 4500], [max(a)*1.2, max(a)*1.2],'color','k')
%plot white flash timing
line([4700 5300], [max(a)*1.2, max(a)*1.2],'color','k')
line([5300 5300], [max(a)*1.2, max(a)*1.3],'color','k')
line([5300 7300], [max(a)*1.3, max(a)*1.3],'color','k')
line([7300 7300], [max(a)*1.2, max(a)*1.3],'color','k')
line([7300 9300], [max(a)*1.2, max(a)*1.2],'color','k')
%fill gap
for i=1:15
    line([4500 4700], [max(a)*1.4/15*(i-1),max(a)*1.4/15*(i)],'color','k')
end
axis([0 9200 0 max(max(a)*1.35,1)])
title('ND4 putative OFF')