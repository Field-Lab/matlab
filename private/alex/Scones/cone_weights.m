
piece = '2012-09-13-2';
run = 'data009';

load('/Volumes/Analysis/2012-09-13-2/subunits/data009_sbc/conepreprocess.mat')

ncells=length(datarun.cell_ids);


% incl - only care for frames with reference cone at desired value
% excl - other cones should be above THR (OFF cells want negative  input)
% only - other cones should EACH have ABS value below THR
% str - other cones SUM should have ABS value below THR
% rand - assign random values, for control

mode='incl'; % 'incl', 'excl', 'rand', 'str', 'onl'
bintotal=linspace(-.5, 0.5,12);
bintotal=bintotal(1:end-1);
ninclude=0; % n of trials to include
plotit=false;
plotbybin=false;
plotalt=true;
plot_each_cone=false;
thr=0.1; % for "rest" in str, only, excl
bkgr=false;
alt_all=cell(1,ncells);
alt_std_all=alt_all;


binwidth=diff(bintotal(1:2));
nbins=length(bintotal);
% alt_coneID=cell(1, length(ncells));
for i = 1:ncells
        
    ncones=size(datarun.cones.weights,1);
%     ncones=sum(datarun.cones.weights(:,i));
%     [val, ic]=sort(conerun.cones.weights(:,cellID), 'descend');
%     alt_coneID{i}=ic(1:ncones(i));
    
%     myinputs=datarun.cone_inputs(:,datarun.cones.weights(:,i)==1);
    myinputs=datarun.cone_inputs;
    myrate=double(datarun.spike_rate(i,:));
    if bkgr
        myrate=myrate-mean(myrate);
    end    

    mycone=zeros(nbins, 26, ncones);
    myconeSTD=mycone;
    alt=zeros(size(mycone,1), ncones);
    alt_std=alt;
    nentries=zeros(ncones,size(mycone,1));
    exclCones=cell(ncones, nbins);
    
    if plotit
        figure
        set(gcf, 'Name', [dates(ncells(i), :), '  Cell ', int2str(rgcs(i)), '    ', mode])
        set(gcf, 'position', [1           1        1920        1105]);
    end
    
    for k=1:ncones
        refcone=k;
        excones=setdiff(1:ncones,k);
        cnt=1;
        for bin=bintotal
            
            switch mode
                case 'incl'
                    tmp=find(myinputs(:,refcone)>= bin & myinputs(:,refcone)< (bin + binwidth));
                case 'excl'
                    tmp=find(myinputs(:,refcone)>= bin & myinputs(:,refcone)< (bin + binwidth) & ...
                        ~sum(myinputs(:,excones)<thr, 2));
                case 'str'
                    rest=sum(myinputs(:,excones),2);
                    tmp=find(myinputs(:,refcone)>= bin & myinputs(:,refcone)< (bin + binwidth) & ...
                        abs(rest)<abs(thr));
                case 'only'
                    tmp=find(myinputs(:,refcone)>= bin & myinputs(:,refcone)< (bin + binwidth) & ...
                        ~sum(abs(myinputs(:,excones))>thr, 2));
                case 'rand'
                    tmp=ceil(rand(100,1)*length(myrate-100));
            end
            
            
            tt=[];
            if ninclude
                trials=ninclude;
            else
                trials=length(tmp);
            end
            if ~isempty(tmp)
                for j=1:trials
                    if tmp(j)>6 && tmp(j)<length(myrate)-21
                        tt=[tt; myrate(tmp(j)-5:tmp(j)+20)];
                    end
                end
                
                
            end
            nentries(k,cnt)=trials;
            mycone(cnt,:,k)=mean(tt);
            %             alt_std(cnt,k)=std(tt(:,7));
            %             mycone(cnt,:,k)=sum(tt);
            alt(cnt,k)=mycone(cnt,7, k);
            
            
            cnt=cnt+1;
        end
        
        if plotit
            subplot(3,ceil(ncones/3),k)
            plot(mycone(:,:,k)')
            title(int2str(nentries(k,:)))
            line([6 6], [-0.1 1.5], 'color', 'k')
            axis tight
            
        end
    end
    
    if plotbybin
        
        figure
        set(gcf, 'Name', [dates(ncells(i), :), '  Cell ', int2str(rgcs(i)) , ' BY BIN   ', mode])
        set(gcf, 'position', [1           1        1920        1105]);
        for j=1:10
            subplot(2,5, j)
            wr=squeeze(mycone(j,:,:));
            plot(wr)
            axis tight
            title(num2str(mean(nentries(:,j))))
        end
    end
    
    if plot_each_cone
        bin=bintotal;
        
        for k=1:ncones(i)
            figure
            set(gcf, 'Name', [dates(ncells(i), :), ', CONE ', int2str(k), ',  Cell ', int2str(rgcs(i)) , ' BY BIN   ', mode])
            set(gcf, 'position', [665   364   644   471]);
            
            mmax=max(mycone(:));
            mmin=min(mycone(:));
            for j=1:nbins
                subplot(5,5, j)
                wr=squeeze(mycone(j,:,k));
                plot(wr, 'r', 'linewidth', 3)
                axis([1 26 mmin mmax])
                line([1 26], [0 0], 'color', 'k')
                title(['[',num2str(bin(j)), ', ', num2str(bin(j)+0.1), '], n=', num2str(nentries(k,j))])
            end
        end
        
        tmp=zeros(ncones(i),26);
        for k=1:ncones(i)
            wr=squeeze(mycone(1:4,:,k));
            
            tmp(k,:)=sum(wr)/sum(nentries(k,1:4));
        end
        figure
        hold on
        plot(tmp', 'linewidth',2)
        title('summed, summed,  mean')
    end
    
    if plotalt
        
        figure
        set(gcf,'Name',['Cell ', int2str(datarun.cell_ids(i)) , ' ALT   ', mode, '  ', int2str(ncones), ' cones'])
        plot(alt, '-x')
%         line([1 10], [0,0],'color','k')
        bin=bintotal;
        set(gca,'xtick',1:nbins, 'xticklabel', {num2str(bin')})
        xlabel('beg of the bin')
        axis tight
        a=get(gca,'YLim');
        line([nbins/2+0.5, nbins/2+0.5], a,'color','k')
        drawnow
    end
    
    alt_all{i}=alt;
    alt_std_all{i}=alt_std;
    
end


figure
thresh=0.25;

for i=1:ncells
    mycrs=alt_all{i}(end-4:end,:);
    mycrs=mycrs(end:-1:1,:)';
    ncones=size(mycrs,1);
    mycrsx=repmat([-0.4:0.1:0],ncones,1);
    mycrs(:,5)=max(mean(mycrs(:,5)), 0);
    
    [p resnorm residual] = normcdfxscalesimple(mycrs, mycrsx, 'plot', false, 'title', false);
    myrfweights = p(1:end-2)';
    myrfweights=myrfweights/max(myrfweights);
    figure
    bar(sort(myrfweights))

end
