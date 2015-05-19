

for i=4:13
    if i<10
        j = ['data00',int2str(i),'-from-d03_13'];
    else 
        j = ['data0',int2str(i),'-from-d03_13'];
    end
    datarun = load_data(['/Volumes/Analysis/2008-03-25-4/d03-13-norefit/',j,'/',j]);
    datarun = load_params(datarun,'verbose',1);
    datarun = load_sta(datarun, 'keep_java_sta', true);
    datarun = set_polarities(datarun);
    datarun = load_neurons(datarun);
    datarun = get_sta_summaries(datarun, 'all');
    my_runs{i-3} =datarun; 
end


filepath = '/Users/alexth/Desktop/Light_adaptation/2008-03-25-4/d03-13-norefit/';
ndfs = [0, 0.6, 1, 1.6, 3.3, 0, 1.3, 2.6, 1.3, 0];
for datarunID = 1:length(datarun.cell_ids)
    visionID = datarun.cell_ids(datarunID);
    [folder, my_type] = find_cell_type(datarun, visionID);
    
    if ~strcmp(folder,'crap') && ~strcmp(folder,'duplicates')
        if ~exist([filepath,folder],'dir')
            mkdir([filepath,folder]);
        end
        
        fig=figure('PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
        set(fig,'color','white','position',[1 24 1886 1074]);
        
        cnt = 1;
        for i=[1 7 6 9 10 4 2 8 3 5]
            datarun = my_runs{i};
            tc = datarun.stas.time_courses{datarunID};
            if datarun.stas.polarities{datarunID}<0
                tmp = -tc;
            else
                tmp = tc;
            end
            tt = max(tmp(:));
            if ~isempty(tt)
                [a,b] = find(tmp == tt, 1); % a is time point, b is color
                
                sta = datarun.stas.stas{datarunID}(:,:,:,a);
                
                tt = 0.5/max(abs(sta(:)));
                
                sta = sta*tt+0.5;
                
                subplot(5,4,cnt)
                imagesc(sta)
                set(gca,'dataaspectratio',[1 1 1])
                
                title(['data00', int2str(i+3), ', NDF', num2str(ndfs(i))])
                
                subplot(5,4,cnt+1)
                
                tc_bin=datarun.stimulus.refresh_period;
                my_bins=-tc_bin*27:tc_bin:tc_bin*2;
                plot(my_bins, tc(:,1), 'r')
                hold on
                plot(my_bins, tc(:,2), 'g')
                plot(my_bins, tc(:,3), 'b')
                line([-1000,200], [0,0],'color','k')
                axis([-1000,100,-Inf,Inf])
                
                if cnt==1
                    title(['data00', int2str(i+3), ', NDF', num2str(ndfs(i)),', cell ', int2str(visionID), ', ',folder]);
                else
                    title(['data00', int2str(i+3), ', NDF', num2str(ndfs(i))])
                end
            end
            cnt = cnt+2;
        end
        
        print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(visionID)]));
        close(fig)
    end
end

