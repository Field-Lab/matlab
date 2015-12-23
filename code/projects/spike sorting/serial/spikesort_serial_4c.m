
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%to do:
%test diffrent spike thresholds
%one template per electrode ... or in neighorhood
%baysian fit
%more robust sig electrode selection
%test whether template is good
%shift fit sometimes broken
%amplitude fit?
%automated duplicate joining

if 0
    clear
    visvolatile();
end

if 1
    if 1
        k=1; 
        
        if 0
            data_path='/Volumes/Twist/Data/Greschner/2011-08-04-4/data000/';
            org_raw_data_path='/Volumes/Twist/Data/Greschner/2011-08-04-4/data000/';
        end
        if 0
            data_path='/Volumes/Twist/Data/Greschner/2011-08-04-4/data000/data000/';
            org_raw_data_path='/Volumes/Twist/Data/Greschner/2011-08-04-4/data000/data000/';
        end
        if 0
            data_path='/Volumes/Stream-phoenix/Data/Greschner/2011-08-04-4/data000';
            org_raw_data_path='/Volumes/Stream-phoenix/Data/Greschner/2011-08-04-4/data000';
        end
        if 0
            data_path='/Volumes/Stream-phoenix/Data/Greschner/2011-08-01-2/data003';
            org_raw_data_path='/Volumes/Stream-phoenix/Data/Greschner/2011-08-01-2/data003';
        end
        if 0
            data_path='/Volumes/Stream-phoenix/MacintoshHD/Data/Greschner/2011-08-04-4/data001';
            org_raw_data_path='/Volumes/Stream-phoenix/MacintoshHD/Data/Greschner/2011-08-04-4/data001';
        end
        if 0
            data_path='/Volumes/Stream-phoenix/MacintoshHD/Data/Greschner/2011-08-01-2/data003';
            org_raw_data_path='/Volumes/Stream-phoenix/MacintoshHD/Data/Greschner/2011-08-01-2/data003';
        end
        if 1
            data_path='/Volumes/Stream-phoenix/Data/Greschner/2011-08-04-1/data000';
            org_raw_data_path='/Volumes/Stream-phoenix/Data/Greschner/2011-08-04-1/data000';
        end

        
        xml_data_path{1}='/snle/lab/Development/RRS/xml-library/xml-library/primate.xml';
        xml_data_path{2}='/snle/lab/Development/RRS/xml-library/xml-library/primate_thr3_5.xml';
        xml_data_path{3}='/snle/lab/Development/RRS/xml-library/xml-library/primate_thr3_0.xml';
        %xml_data_path{4}='/snle/lab/Development/RRS/xml-library/xml-library/primate_thr2_5.xml';
        xml_switch=[0 25 50];
 
        
        array=512;%519;
        elec_nr=array+1;
        
        
        samplingrate=20000;
        
        raw_data_path=org_raw_data_path;
        
        debug=1;
        
    end
    
    for k=k:99
        output_data_path=sprintf('%s/serial_%02d/',data_path,k);
    
        
        display(sprintf('\n\n%% serial_%02d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',k));
        display(sprintf('\n%s\n\ninput: %s\noutput: %s\n',datestr(now),raw_data_path,output_data_path))
        tic
        
        
        if 1 %run vision
            if ~exist(output_data_path,'dir')
                command=sprintf('/snle/lab/Development/scripts/vision-auto-neuronidentification %s %s false true %s',raw_data_path,output_data_path,xml_data_path{max(find(xml_switch<k))});
                unix(command);
            end
        end
    
       
        
        if 1 %load files
            clear datarun
            t=dir([output_data_path '/*.neurons-raw']);
                datarun.names.rrs_neurons_path=[output_data_path '/' t.name];
            t=dir([output_data_path '/*.globals']);
                datarun.names.rrs_globals_path=[output_data_path '/' t.name]; 
            t=dir([output_data_path '/*.prj']);
                datarun.names.rrs_projection_path=[output_data_path '/' t.name]; 
            t=dir([output_data_path '/*.model']);
                datarun.names.rrs_model_path=[output_data_path '/' t.name]; 

            opt=struct('verbose',1,'load_neurons',1);
            datarun=load_data(datarun,opt);    

            prjObject = edu.ucsc.neurobiology.vision.matlab.ReadProjections(datarun.names.rrs_projection_path);
            model = edu.ucsc.neurobiology.vision.io.ClusteringModelFile(datarun.names.rrs_model_path);
        end
    
        
    
        if 1 %find good cluster
            thr_contam=0;

            clear tcluster
            a=1;
            T=text_waitbar('find clusters');
            for j=1:elec_nr-1
                
                prjObject.readElectrode(j);
                prj=prjObject.getProjections();
                nSpikes=prjObject.getSpikeCount();
                prjs=prj(:,1:nSpikes);

                clusters = model.getNeuronExtraction(j);
                
                if ~isempty(clusters)
                    like=[];
                    for i=1:size(clusters.means,1)
                        if any(clusters.covariances(i,:))
                            gg=gmdistribution(clusters.means(i,:),clusters.covariances(i,:));
                            like(:,i)=pdf(gg,prjs');%*clusters.probability(i);
                        end
                    end

                    if 0 %plot   
                        clf 
                        di=[1 2;1 3;2 3];
                        col=colormap('lines');
                        [junk tlike]=max(like,[],2);
                        for i=1:size(clusters.means,1)
                            for ii=1:3
                                subplot(1,3,ii)
                                hold on
                                tt=find(tlike==i);
                                tt=tt(1:500);
                                plot(prjs(di(ii,1),tt),prjs(di(ii,2),tt),'.','color',col(i,:));
                            end
                        end
                    end
                    
                    if size(like,2)>1
                        [s1 s2]=sort(like,2);
                        t_like=s2(:,end);
                        tr_like=s1(:,end)./s1(:,end-1);

                        for i=1:size(clusters.means,1)
                            t1=find(t_like==i);

                            cell_id=(j-1)*15+i;

                            index=get_cell_indices(datarun,cell_id);
                            if ~isempty(datarun.spikes{index})
                                datarun=get_contamination(datarun,cell_id);

                                if(datarun.contamination(index)<=thr_contam);
                                    tcluster(a).cell_id=cell_id;
                                    tcluster(a).electrode=j;
                                    tcluster(a).cluster=i;
                                    tcluster(a).r_like=robust_mean(log(tr_like(t1)));
                                    a=a+1;
                                end
                            end
                        end
                    end
                end
                T=text_waitbar(T,j/(elec_nr-1));
            end
            
            if 0
                h=[tcluster.r_like];
                hist(h(find(isfinite(h))),100);
            end

        end
        
  
       
        if 1 %select templates 
            dup_thr=.3;
            nr_neighbor=1;
            
            display('select clusters...');
            
            templates=[];
            
            cell_ids=[tcluster.cell_id];
            r_like=[tcluster.r_like];
  
            display(sprintf('cluster found: %d',length(cell_ids)))
   
            %greedy cluster selection
            [junk sortind]=sort(r_like,'descend');
            
                %one per neighborhood
                if 0
                    m=zeros(elec_nr,1);
                    a=1;
                    for i=1:length(sortind)

                        %one per neighborhood
                        tneighbors=get_ei_neighbors(tcluster(sortind(i)).electrode,array,nr_neighbor);
                        if all(~m(tneighbors))

                            %check duplicate
                            if a>1
                                tcell_ids=[tcluster(sortind(i)).cell_id templates.template.cell_id];
                                datarun=correlation(datarun, tcell_ids);
                                tind=get_cell_indices(datarun,tcell_ids);
                                syn=datarun.synchrony_index(tind(1),tind(2:end));
                            else
                                syn=0;
                            end

                            if syn<dup_thr
                                m(tneighbors)=1;
                                templates.template(a)=tcluster(sortind(i));
                                a=a+1;
                            end
                        end
                    end
                end
            
                %best 100
                if 0
                    nr_cluster=75;

                    a=1;
                    i=1;
                    dup=0;
                    while a<=nr_cluster & i<=length(sortind)
                        %check duplicate
                        if a>1
                            tcell_ids=[tcluster(sortind(i)).cell_id templates.template.cell_id];
                            datarun=correlation(datarun, tcell_ids);
                            tind=get_cell_indices(datarun,tcell_ids);
                            syn=datarun.synchrony_index(tind(1),tind(2:end));
                        else
                            syn=0;
                        end

                        if all(syn<dup_thr)
                            templates.template(a)=tcluster(sortind(i));
                            a=a+1;
                        else
                            dup=dup+1;
                        end
                        
                        i=i+1;
                    end
                end               
                
                
                %best 100, same elec
                if 1
                    nr_cluster=40;

                    a=1;
                    i=1;
                    m=zeros(elec_nr,1);
                    while a<=nr_cluster & i<=length(sortind)
                        %check duplicate
                        if a>1
                            tcell_ids=[tcluster(sortind(i)).cell_id templates.template.cell_id];
                            datarun=correlation(datarun, tcell_ids);
                            tind=get_cell_indices(datarun,tcell_ids);
                            syn=datarun.synchrony_index(tind(1),tind(2:end));
                        else
                            syn=0;
                        end

                        if all(syn<dup_thr) & ~m(tcluster(sortind(i)).electrode)
                            templates.template(a)=tcluster(sortind(i));
                            a=a+1;
                            m(tcluster(sortind(i)).electrode)=1;
                        end
                        
                        i=i+1;
                    end
                end    
                
                
                
                %std of likelihood
                if 0
                    std_like_thr=5;
                    
                    h=[tcluster.r_like];
                    h(find(~isfinite(h)))=max(isfinite(h));
                    bin=[min(h):(max(h)-min(h))/200:max(h)];
                    hh=histc(h,bin);
                    [s,m,a]=gaussfit(bin,hh,.2);
                    tstd_like_thr=m+s*std_like_thr;
                    
                    if 1%plot
                        clf 
                        subplot(2,1,1)
                            hold on
                            plot(bin,hh)
                            y=a*exp( -(bin-m).^2 / (2*s^2) );
                            plot(bin,y,'r')
                            plot(tstd_like_thr,1,'r*')
                        subplot(2,1,2)
                            hold on
                            plot(sort(h,'descend'),'.-')
                            plot([0 500],[1 1]*tstd_like_thr,'r')
                            
                        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 20 10]); 
                        print('-djpeg',sprintf('%s/likelihood.jpg',output_data_path), '-r100');  
                    end

                    
                    a=1;
                    i=1;
                    dup=0;
                    m=zeros(elec_nr,1);
                    while tcluster(sortind(i)).r_like>tstd_like_thr & i<=length(sortind)
                        %check duplicate
                        if a>1
                            tcell_ids=[tcluster(sortind(i)).cell_id templates.template.cell_id];
                            datarun=correlation(datarun, tcell_ids);
                            tind=get_cell_indices(datarun,tcell_ids);
                            syn=datarun.synchrony_index(tind(1),tind(2:end));
                        else
                            syn=0;
                        end

                        if all(syn<dup_thr) & ~m(tcluster(sortind(i)).electrode)
                            templates.template(a)=tcluster(sortind(i));
                            a=a+1;
                            m(tcluster(sortind(i)).electrode)=1;
                        else
                            dup=dup+1;
                        end
                        
                        i=i+1;
                    end
                end   
                
                
                
            if 0 %plot
                clf
                hold on
                h=[tcluster.r_like];
                hist(h(find(isfinite(h))),300);
                h = findobj(gca,'Type','patch');
                set(h,'FaceColor','r','EdgeColor','w')
                
                h=[templates.template.r_like];
                hist(h(find(isfinite(h))),300);
            end
            
            templates.cell_ids=[templates.template.cell_id];
            
            display(sprintf('cluster selected: %d',length(templates.cell_ids))) 
            display(sprintf('median cluster likelihood: %.2f (min:%.2f max:%.2f)\n',median([templates.template.r_like]),min([templates.template.r_like]),max([templates.template.r_like])))

            index=get_cell_indices(datarun,templates.cell_ids);
            templates.sptimes=[];
            templates.template_index=[];  
            for i=1:length(templates.cell_ids)
                templates.sptimes=[templates.sptimes; datarun.spikes{index(i)}];
                templates.template_index=[templates.template_index; ones(size(datarun.spikes{index(i)}))*i]; 
            end
  
        end

       
        
        if 1 %fast ei calculation - make templates use ei
            nlPoints=30;
            nrPoints=90;
            spike_nr_template=5000;
            sig_thr1=30;%pick elec for templates
            sig_thr2=0.5;%pick elec for fit
            sig_thr3=2;%mask in time
 
            shift_bin=.1;
            tshift=[-6:6];
            
            if ~exist([output_data_path '/temp/temp.ei'])
                display('calculate EIs')
                % save temp.neuron
                    old_neurons = edu.ucsc.neurobiology.vision.io.NeuronFile(datarun.names.rrs_neurons_path);
                    header = old_neurons.getHeader;
                    triggers = old_neurons.getTTLTimes;
                    old_neurons.close;

                    index=get_cell_indices(datarun,templates.cell_ids);
                    clear spikes
                    for i=1:length(index)
                        spikes{i}=int32(datarun.spikes{index(i)}*samplingrate); 
                    end

                    unix(['mkdir ' output_data_path '/temp']);
                    unix(['cp ' output_data_path '/*.globals ' output_data_path '/temp/temp.globals']);

                    save_neurons([output_data_path '/temp/temp.neurons'], header, triggers, spikes, datarun.cell_ids(index), datarun.channels(index));

                    command=sprintf('/snle/lab/Development/scripts/vision-auto-ei-only-64 %s %s/temp %d %d %d',raw_data_path,output_data_path, nlPoints, nrPoints, spike_nr_template);
                    unix(command);
            end
            
            %load ei
                clear dataset
                dataset.names.rrs_ei_path=[output_data_path '/temp/temp.ei'];
                dataset.names.rrs_globals_path=[output_data_path '/temp/temp.globals'];
                dataset.names.rrs_neurons_path=[output_data_path '/temp/temp.neurons'];
                opt=struct('verbose',1,'load_neurons',1,'load_ei',1);
                dataset=load_data(dataset,opt);    

            T=text_waitbar('make templates');
            for j=1:length(templates.cell_ids)
                
                index=get_cell_indices(dataset,templates.cell_ids(j));
                ei=dataset.ei.eis{index};
                ei(:,2:end)=ei(:,1:end-1);%hack match old code
                
                %thresholds
                [t1]=(max(abs(ei),[],2));
                [junk t2]=sort(t1);
                t3=ei(t2(1:end-5),1:20);
                tstd=robust_std(t3(:));
                sig_elec=find(t1>tstd*sig_thr1);
                if length(sig_elec)<2
                    sig_elec=t2(end-1:end);
                end
                
                [ttt1 ttt2]=sort(t1);
                if ttt1(end-1)>ttt1(end)*sig_thr2
                    sig_elec_fit=ttt2(end-1:end);
                else
                    sig_elec_fit=ttt2(end);
                end 
                

                
                ttemplate=ei(sig_elec,:);
                
                %mask in time
                for i=1:length(sig_elec)
                    
                    %does not work reliable ttemplate(i,:)=ttemplate(i,:)-mean(ttemplate(i,1:15));
                    
                    t=circshift(conv(abs(ttemplate(i,:)),ones(1,20)/20,'same')',-10);
                    tt=find(t<tstd*sig_thr3);
                    ttemplate(i,tt)=0;
                    
                    if 0%plot
                        clf
                        hold on
                        plot(ttemplate(i,:)')
                        plot(t','g')
                        rr(i,tt)=0;
                        plot(rr(i,:)','r')
                        pause
                    end 
                end
                
                templates.template(j).sig_elec=sig_elec+1;%shift for trigger channel
                templates.template(j).sig_elec_fit=sig_elec_fit+1;
                for i=1:length(templates.template(j).sig_elec_fit)
                    templates.template(j).sig_elec_fit(i,2)=find(templates.template(j).sig_elec==templates.template(j).sig_elec_fit(i,1));
                end
                
                %shift templates
                    x=[1:size(ttemplate,2)];
                    xi=[1:shift_bin:size(ttemplate,2)]; 
                    ti=interp1(x,ttemplate',xi,'cubic');
                    for i=1:length(tshift)
                        tti=circshift(ti,tshift(i));
                        ttii=interp1(xi,tti,x,'cubic');
                        templates.template(j).shift(i).mean=int16(ttii(1:end-1,:));
                        if 0
                            hold on
                            plot(templates.template(j).shift(i).mean)
                        end
                    end
                
                    
                if 0 %plot ei
                    clf;
                    subplot(2,1,1)
                        positions = electrode_positions(array);
                        plot_ei_(ei,positions,0,'scale',5,'flipud',0,'cutoff',0.0,'sub_scale',.0) 
                        hold on
                        plot(positions(sig_elec,1),positions(sig_elec,2),'r.')
                        plot(positions(templates.template(j).electrode,1),positions(templates.template(j).electrode,2),'r*')
                        plot(positions(sig_elec_fit,1),positions(sig_elec_fit,2),'g.')
                    subplot(2,1,2)
                        plot(templates.template(j).shift(1).mean)
                        
                    pause
                end 

                if debug %plots
                    clf;
                    subplot(1,2,1)
                        positions = electrode_positions(array);
                        plot_ei_(ei,positions,0,'scale',5,'flipud',0,'cutoff',0.0,'sub_scale',.0) 
                        hold on
                        plot(positions(sig_elec,1),positions(sig_elec,2),'r.')
                        plot(positions(templates.template(j).electrode,1),positions(templates.template(j).electrode,2),'r*')
                        plot(positions(sig_elec_fit,1),positions(sig_elec_fit,2),'g.')
                    subplot(1,2,2)
                        plot(templates.template(j).shift(1).mean)
                    %saveas(gcf,sprintf('%s/template_%03d_%03d_%04d.jpg',output_data_path,j,templates.template(j).electrode,templates.template(j).cell_id),'jpg')   
                    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 20 10]); 
                    print('-djpeg',sprintf('%s/template_%03d_%03d_%04d.jpg',output_data_path,j,templates.template(j).electrode,templates.template(j).cell_id), '-r100');
                end
                
                T=text_waitbar(T,j/length(templates.cell_ids));
            end
            
            t=zeros(size(templates.cell_ids));
            for i=1:length(templates.cell_ids)
                t(i)=length(templates.template(i).sig_elec);
            end 
            
            display(sprintf('median significant electrodes: %d (min:%d max:%d)',median(t),min(t),max(t)))   

        end
        
 
        
 
        if 1 %save 
            display(['save ' output_data_path '/templates.mat']);
            save([output_data_path '/templates.mat'],'templates');
        end
        
        %test [s1 s2]=sort(templates.sptimes);templates.sptimes=templates.sptimes(s2);templates.template_index=templates.template_index(s2);
        
        
        
        
        if 1 %substract in raw data
            verbose=0;
            
            %re-use old header
            if ~exist('rdfh')
                t=dir([org_raw_data_path '/*.bin']);              
                rdf=edu.ucsc.neurobiology.vision.io.RawDataFile([org_raw_data_path '/' t(1).name]);
                rdfh=rdf.getHeader();
            end

            rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(raw_data_path);
            
            rds=rdsopen([output_data_path '/bin/'], rdfh);

            T=text_waitbar('process bin');
            tstep=10*samplingrate;
            tborder=int16(zeros(nrPoints,elec_nr));
            tsptimes=round(templates.sptimes*samplingrate);
            h=zeros(size(tshift));
            for i=1:datarun.duration*samplingrate/tstep;

                
                m=rawFile.getData((i-1)*tstep, tstep);
                
                
                %substract leading border
                    m(1:nrPoints,:)=m(1:nrPoints,:)-tborder; %substract leftover from last data piece
                    
                    t=find(tsptimes>=(i-1)*tstep & tsptimes<(i-1)*tstep+nlPoints);
                    for ii=1:length(t)
                        tt=tsptimes(t(ii))-(i-1)*tstep;
                        ttem=templates.template_index(t(ii));
                        
                        if verbose %plot
                            clf
                            for iii=1:min([length(templates.template(ttem).sig_elec) 16])
                                subplot(4,4,iii)
                                hold on
                                plot(m(1:tt+nrPoints,templates.template(ttem).sig_elec(iii)));
                            end
                        end
                        
                        for iii=1:length(tshift)
                            h1=m(1:tt+nrPoints,templates.template(ttem).sig_elec_fit(:,1));
                            h2=templates.template(ttem).shift(iii).mean(end-size(h1,1)+1:end,templates.template(ttem).sig_elec_fit(:,2));
                            h(iii)=rmse(h1(:),h2(:));
                        end
                        [junk tsh]=min(h);

                        m(1:tt+nrPoints,templates.template(ttem).sig_elec)=m(1:tt+nrPoints,templates.template(ttem).sig_elec)-templates.template(ttem).shift(tsh).mean(end-size(h1,1)+1:end,:);
                        
                        if verbose %plot
                            for iii=1:min([length(templates.template(ttem).sig_elec) 16])
                                subplot(4,4,iii)
                                plot(m(1:tt+nrPoints,templates.template(ttem).sig_elec(iii)),'r');
                            end
                            subplot(4,4,1)
                            title('leading border')
                            pause
                        end
                    end 
                    
                     
                %substract trailing border %ignor leftover from following data piece
                    t=find(tsptimes>=i*tstep-nrPoints & tsptimes<i*tstep); 
                    for ii=1:length(t)
                        tt=tsptimes(t(ii))-(i-1)*tstep;
                        ttem=templates.template_index(t(ii));
                        
                        if verbose %plot
                            clf
                            for iii=1:min([length(templates.template(ttem).sig_elec) 16])
                                subplot(4,4,iii)
                                hold on
                                plot(m(tt-nlPoints+1:end,templates.template(ttem).sig_elec(iii)));
                                plot(templates.template(ttem).shift(10).mean(:,iii),'g');
                            end
                        end
                        
                        for iii=1:length(tshift) %make fit only here, even if spike actually in next data piece
                            h1=m(tt-nlPoints+1:end,templates.template(ttem).sig_elec_fit(:,1));
                            h2=templates.template(ttem).shift(iii).mean(1:size(h1,1),templates.template(ttem).sig_elec_fit(:,2));
                            h(iii)=rmse(h1(:),h2(:));
                        end
                        [junk tsh]=min(h);

                        m(tt-nlPoints+1:end,templates.template(ttem).sig_elec)=m(tt-nlPoints+1:end,templates.template(ttem).sig_elec)-templates.template(ttem).shift(tsh).mean(1:size(h1,1),:);
                        t3=templates.template(ttem).shift(tsh).mean(size(h1,1)+1:end,:);
                        tborder(1:size(t3,1),templates.template(ttem).sig_elec)=tborder(1:size(t3,1),templates.template(ttem).sig_elec)+t3;
                        
                        if verbose %plot
                            for iii=1:min([length(templates.template(ttem).sig_elec) 16])
                                subplot(4,4,iii)
                                plot(m(tt-nlPoints+1:end,templates.template(ttem).sig_elec(iii)),'r');
                                %plot(templates.template(ttem).shift(10).mean(:,iii),'g');
                            end
                            subplot(4,4,1)
                            title('trailing border')
                            pause
                        end
                    end
                    
                
                    
                %substract in raw data
                    t=find(tsptimes>=(i-1)*tstep+nlPoints & tsptimes<i*tstep-nrPoints);
                    
                    for ii=1:length(t)
                        tt=tsptimes(t(ii))-(i-1)*tstep;
                        ttem=templates.template_index(t(ii));
                        
                        if verbose %plot
                            clf
                            for iii=1:min([length(templates.template(ttem).sig_elec) 15])
                                subplot(4,4,iii)
                                hold on
                                plot(m(tt-nlPoints+1:tt+nrPoints,templates.template(ttem).sig_elec(iii)));
                            end
                        end
                        
                        for iii=1:length(tshift)
                            h1=m(tt-nlPoints+1:tt+nrPoints,templates.template(ttem).sig_elec_fit(:,1));
                            h2=templates.template(ttem).shift(iii).mean(:,templates.template(ttem).sig_elec_fit(:,2));
                            h(iii)=rmse(h1(:),h2(:));
                        end
                        [junk tsh]=min(h);

                        m(tt-nlPoints+1:tt+nrPoints,templates.template(ttem).sig_elec)=m(tt-nlPoints+1:tt+nrPoints,templates.template(ttem).sig_elec)-templates.template(ttem).shift(tsh).mean;
                        
                        if verbose %plot
                            for iii=1:min([length(templates.template(ttem).sig_elec) 15])
                                subplot(4,4,iii)
                                plot(m(tt-nlPoints+1:tt+nrPoints,templates.template(ttem).sig_elec(iii)),'r');
                            end
                            subplot(4,4,16)
                            plot(h)
                            hold on
                            plot(tsh,h(tsh),'r.')
                            subplot(4,4,1)
                            title(sprintf('%d/%d',ii,length(t)))
                            pause
                        end
                    
                    end

                edu.ucsc.neurobiology.vision.matlab.Matlab.saveRawData(rds, m);
                
                T=text_waitbar(T,i/(datarun.duration*samplingrate/tstep));
            end

            rdsclose(rds);            
            rawFile.close();
        end

        
        if 1 %clean up and start over
            unix(['mkdir ' output_data_path '/bin/data999']);
            unix(['mv ' output_data_path '/bin/*/*.bin ' output_data_path '/bin/data999/data999000.bin']);
            raw_data_path=[output_data_path '/bin/data999/'];
            
            del_files{1}='cov'; del_files{2}='ncov'; del_files{3}='noise'; del_files{4}='params'; del_files{5}='spikes'; del_files{6}='wcov';
            for i=1:length(del_files)
                unix(['rm ' output_data_path '/*.' del_files{i}]);
            end
            
            try
                unix(sprintf('rm -r %s/serial_%02d/bin/',data_path,k-2));
            catch
            end
                
            %unix(['rm -r ' output_data_path '/temp']);
        end
        
        toc
    end

end












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0 %save neurons
    
    %data_path='/Volumes/Twist/Data/Greschner/2011-08-04-4/data000/ver1';
    %data_path='/Volumes/Stream-phoenix/Data/Greschner/2011-08-04-4/data000';
    data_path='/Volumes//Data/Greschner/2011-08-01-2/data003';
    %data_path='/Volumes/Stream-phoenix/MacintoshHD/Data/Greschner/2011-08-04-4/data000';
    
    nr=99;
    
    import('edu.ucsc.neurobiology.vision.matlab.*');
    import('edu.ucsc.neurobiology.vision.io.*');
   
    tnames=dir([data_path '/serial_*']);
    
    sampling_frequency=edu.ucsc.neurobiology.vision.util.VisionParams.samplesPerMillisecond * 1000;
    
    new_neurons_file = [data_path '/serial.neurons'];
    
    %load data
    if 0
        clear res rres;
        cell_ids=[];
        electrodes=[];
        trate=[];
        tlike=[];
        tamp=[];
        tnrcluster=[];
        for k=1:min([length(tnames) nr])-1

            %load files
            clear datarun
            t=dir([data_path '/' tnames(k).name '/*.neurons-raw']);
                datarun.names.rrs_neurons_path=[data_path '/' tnames(k).name '/' t.name];
            t=dir([data_path '/' tnames(k).name '/*.globals']);
                datarun.names.rrs_globals_path=[data_path '/' tnames(k).name '/' t.name]; 

            opt=struct('verbose',1,'load_neurons',1);
            datarun=load_data(datarun,opt);   

            %load temp    
            load([data_path '/' tnames(k).name '/templates.mat']);

            index=get_cell_indices(datarun,templates.cell_ids);
            res(length(cell_ids)+1:length(index)+length(cell_ids))=datarun.spikes(index);
            
            tcell_ids=templates.cell_ids+10000*k;
            cell_ids=[cell_ids tcell_ids];
            rres{k}=tcell_ids;
            electrodes=[electrodes [templates.template.electrode] ];  
            display(length(cell_ids))
            
            for i=1:length(templates.cell_ids)
                tnrcluster=[tnrcluster length(datarun.cell_ids)];
            	trate=[trate length(datarun.spikes{index(i)})];
                tlike=[tlike templates.template(i).r_like];
                tamp=[tamp max(abs(templates.template(i).shift(1).mean(:)))];
            end
        end 
        if 0
            k=min([length(tnames) nr])

            %load files
            clear datarun
            t=dir([data_path '/' tnames(k).name '/*.neurons']);
                datarun.names.rrs_neurons_path=[data_path '/' tnames(k).name '/' t.name];
            t=dir([data_path '/' tnames(k).name '/*.globals']);
                datarun.names.rrs_globals_path=[data_path '/' tnames(k).name '/' t.name]; 

            opt=struct('verbose',1,'load_neurons',1);
            datarun=load_data(datarun,opt);   

            tcell_ids=[length(cell_ids)+1:+length(datarun.spikes)+length(cell_ids)];
            res(tcell_ids)=datarun.spikes;
            cell_ids=[cell_ids tcell_ids];
            rres{k}=tcell_ids;
            electrodes=[electrodes datarun.channels'];   
            display(length(cell_ids))
        end  
    end
    
    %write neuron
    if 0 
        % grab/set old header information
            old_neurons = NeuronFile(datarun.names.rrs_neurons_path);
            header = old_neurons.getHeader;
            triggers = old_neurons.getTTLTimes;

            clear spikes
            for i=1:length(res)
                spikes{i}=int32(res{i}*sampling_frequency); 
            end

        % save out the neurons file
            save_neurons(new_neurons_file, header, triggers, spikes, cell_ids, electrodes);
    end
    
    %write classification.txt
    if 0
        [fid1,msg]=fopen([data_path '/serial-classification.txt'], 'w');
        for i=1:length(rres)
            for ii=1:length(rres{i})
                fprintf(fid1,'%d  All/loop_%02d\n',rres{i}(ii),i);
            end
        end
        fclose(fid1);
    end
    
    if 1
        clf; set(gcf, 'color', 'white');

        ttlike=tlike;
        ttlike(find(~isfinite(ttlike)))=max(ttlike(find(isfinite(ttlike))));
        ttlike=ttlike/max(ttlike);

        subplot(3,1,1)
            hold on
            plot(trate/datarun.duration)
            plot(ttlike*max(trate/datarun.duration),'r');
            ylabel('rate')

        subplot(3,1,2)
            hold on
            plot(tamp)
            plot(ttlike*double(max(tamp)),'r');
            ylabel('amp')
            
        subplot(3,1,3)
            hold on
            t1=zeros(length(rres),1);
            t2=zeros(length(rres),1);n=0;
            for i=1:length(rres)
                t1(i)=length(rres{i});
                t2(i)=robust_mean(tlike(n+1:n+length(rres{i})));
                n=n+length(rres{i});
            end
            plot(t1)
            plot(t2/max(t2)*max(t1),'r')
    end
  
end
%vision-auto-sta-only-64 /Volumes/Stream-phoenix/Data/Greschner/2011-08-04-4/data000/serial /snle/acquisition/movie-xml/RGB-8-16-0.48-11111-80x60.xml






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make cell-type summery plot
if 0
    if 0
        clear datarun
        datarun.names.rrs_classification_path='/Volumes/Twist/Data/Greschner/2011-08-04-4/data000/ver1/serial/serial.classification.txt';
        datarun = load_vision_cell_types(datarun, varargin);
    end
    
    cell_ids=[];
    for i=1:length(datarun.cell_types)
        cell_ids=[cell_ids datarun.cell_types{i}.cell_ids];
    end
    
    r=zeros(length(datarun.cell_types),length(cell_ids));
    for i=1:length(datarun.cell_types)
        if ~isempty(datarun.cell_types{i}.cell_ids)
            t=cumsum(histc(datarun.cell_types{i}.cell_ids,[1:length(cell_ids)]));
            r(i,:)=t/max(t);
        end
    end
    
    clf; set(gcf, 'color', 'white');
    hold on
    n=1:5;
    plot(r(n,:)');
    clear l
    for i=1:length(n)
        l{i}=datarun.cell_types{n(i)}.name;
    end
    legend(l,2);
    box on
    
end









%ei sub test
if 0
    if 0
        display('calculate EIs')

        unix(['mkdir ' output_data_path '/temp/temp']);
        unix(['cp ' output_data_path '/temp/temp.globals ' output_data_path '/temp/temp/temp.globals']);
        unix(['cp ' output_data_path '/temp/temp.neurons ' output_data_path '/temp/temp/temp.neurons']);

        command=sprintf('/snle/lab/Development/scripts/vision-auto-ei-only-64 %s %s/temp/temp %d %d %d',[output_data_path '/bin/data999/'],output_data_path, nlPoints, nrPoints, spike_nr_template);
        unix(command);
    end
    if 0
        %load ei
        clear d1
        d1.names.rrs_ei_path=[output_data_path '/temp/temp.ei'];
        d1.names.rrs_globals_path=[output_data_path '/temp/temp.globals'];
        d1.names.rrs_neurons_path=[output_data_path '/temp/temp.neurons'];
        opt=struct('verbose',1,'load_neurons',1,'load_ei',1);
        d1=load_data(d1,opt); 

        clear d2
        d2.names.rrs_ei_path=[output_data_path '/temp-test/temp/temp.ei'];
        d2.names.rrs_globals_path=[output_data_path '/temp-test/temp/temp.globals'];
        d2.names.rrs_neurons_path=[output_data_path '/temp-test/temp/temp.neurons'];
        opt=struct('verbose',1,'load_neurons',1,'load_ei',1);
        d2=load_data(d2,opt); 
    end
    
    j=8;
    index=get_cell_indices(d1,templates.cell_ids(j));
    ei1=d1.ei.eis{index};
    index=get_cell_indices(d2,templates.cell_ids(j));
    ei2=d2.ei.eis{index};
    
    yr=[min(ei1(:)) max(ei1(:))];
    clf;
    subplot(2,2,1)
    plot(ei1')
    set(gca, 'ylim', yr)
    subplot(2,2,2)
    plot(ei2')
    set(gca, 'ylim', yr)
    subplot(2,2,3)
    plot(ei1(templates.template(j).sig_elec-1,:)')
    set(gca, 'ylim', yr)
    subplot(2,2,4)
    plot(ei2(templates.template(j).sig_elec-1,:)')
    set(gca, 'ylim', yr)

end    






















