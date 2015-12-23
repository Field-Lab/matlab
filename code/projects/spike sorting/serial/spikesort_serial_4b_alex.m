
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

%start messing around on line 400      nlPoints=30; nrPoints=90;spike_nr_template=5000 sig_thr1=30;%pick elec for templates sig_thr2=0.5;%pick elec for fit sig_thr3=2;%mask in time
%line 540 ish is verbose..  herr you can see the affects of template
%subtraction


if 1  %setting parameters of what sort?
    clear
    visvolatile();
end
totaliter=10;
totalcluster=800;
importantinfo=zeros(

if 1
    if 1
        
        %where to look  for the data
        if 1
            data_path='/Volumes/Wheel/Data/Akheitman/PlayingAround/data000';  %string
            org_raw_data_path='/Volumes/Wheel/Data/Akheitman/PlayingAround/data000';%original data path
        end

        
        %Parameters with different voltage thresholds  used in the commmand
        %of actually running the neuroidentification  - needs 
        %vision-auto-neuronidentification <input-dir-full-path> <output-dir-full-path> [ei] [whitened] [config-xml-file]
        %
        xml_data_path{1}='/snle/lab/Development/RRS/xml-library/xml-library/primate.xml';
        xml_data_path{2}='/snle/lab/Development/RRS/xml-library/xml-library/primate_thr3_5.xml';
        xml_data_path{3}='/snle/lab/Development/RRS/xml-library/xml-library/primate_thr3_0.xml';
        xml_switch=[0 5 8];
 
        elec_nr=520;
        array=519;
        samplingrate=20000;
        
        raw_data_path=org_raw_data_path;
        
        debug=1;
        
    end
    
    for k=1:totaliter
        output_data_path=sprintf('%s/serial_%02d/',data_path,k); % SPRINTF Write formatted data to string.
    
        
        display(sprintf('\n\n%% serial_%02d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',k));% displays serial_iterationnumber
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
    
        
    
        if 1 %find good cluster   everythingsalready automatically clustered by algoritm
            thr_contam=0;

            clear tcluster
            a=1;   %counts uncontamiated cell
            T=text_waitbar('find clusters');   %even though clusters are already found in a sense
            for j=1:elec_nr-1
                
                prjObject.readElectrode(j);   %is this an execution...   what does this 1 or 0 do?
                prj=prjObject.getProjections();    % this isn't changing how it should???
                nSpikes=prjObject.getSpikeCount();   % this is changig how it should??
                prjs=prj(:,1:nSpikes);

                clusters = model.getNeuronExtraction(j);
                
                if ~isempty(clusters) 
                    like=[];
                    for i=1:size(clusters.means,1)  %the number of clusters
                        if any(clusters.covariances(i,:))   % any is a boolean that sees if there exists a non zero element.  (true if there exists,  falso if doesn't)
                            gg=gmdistribution(clusters.means(i,:),clusters.covariances(i,:));
                            like(:,i)=pdf(gg,prjs');%*clusters.probability(i);  still going electrode by electrode
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
                    
                    [s1 s2]=sort(like,2);  % s2 is the index matrix that keeps track cluster number
                    t_like=s2(:,end);  %identifies the higheset prob cluster for each spike
                    tr_like=s1(:,end)./s1(:,end-1);   % the "gain " in probability  (the length of number of spikes)   should all be greater than 1  bigger the better
                    
                    for i=1:size(clusters.means,1)  %number of clusters
                        t1=find(t_like==i);  %a series of indices

                        cell_id=(j-1)*15+i;  % index cell id
                        
                        index=get_cell_indices(datarun,cell_id);
                        if ~isempty(datarun.spikes{index})
                            datarun=get_contamination(datarun,cell_id);

                            if(datarun.contamination(index)<=thr_contam); % contamination check  brings cell count down... and we get info for each on tcluster
                                tcluster(a).cell_id=cell_id;
                                tcluster(a).electrode=j;
                                tcluster(a).cluster=i;
                                tcluster(a).r_like=robust_mean(log(tr_like(t1)));        %  robust mean of the log of the "gain"
                                a=a+1;  %counting a non-contam cluster
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
        
  
       
        if 1 %select templates chossing which clusters of the the non-contam we will use for subtraction
            dup_thr=.3;
            nr_neighbor=1;
            
            display('select clusters...');
            
            templates=[];
            
            cell_ids=[tcluster.cell_id];   %vectorizing   tcluster data
            r_like=[tcluster.r_like];   % we should ideally chose the best from these
  
            display(sprintf('cluster found: %d',length(cell_ids)))   %showing how many non-contaminted clusters we found
   
            %greedy cluster selection
            [junk sortind]=sort(r_like,'descend');   % junk are the r values from best to worst... sortind are the cluster index from best to worst
            
 
                
            if 1   % now looking to chose 40 best clusters from which to undertake subtraction
               nr_cluster=40;  %don't want tooo large

               a=1;
               i=1;
               m=zeros(elec_nr,1);
               while a<=nr_cluster & i<=length(sortind) 
                   %check duplicate
                   if a>1
                       tcell_ids=[tcluster(sortind(i)).cell_id templates.template.cell_id];  %cell ids always greater than or eq to cluster index
                       datarun=correlation(datarun, tcell_ids);
                       tind=get_cell_indices(datarun,tcell_ids);  
                       syn=datarun.synchrony_index(tind(1),tind(2:end));   %datarun goes by cluster index so need tind
                   else
                       syn=0;
                   end
                   
                   if all(syn<dup_thr) & ~m(tcluster(sortind(i)).electrode)  % we are not being duplicated by any of the other top 40 and we are not on the same electrode
                       templates.template(a)=tcluster(sortind(i));
                       a=a+1;
                       m(tcluster(sortind(i)).electrode)=1;  %m=1  means the electrode is important,makrs it a 1 instead of 0   
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
            
            templates.cell_ids=[templates.template.cell_id];   %somehwo one shifted from tcell_ids
            
            display(sprintf('cluster selected: %d',length(templates.cell_ids))) 
            display(sprintf('median cluster likelihood: %.2f (min:%.2f max:%.2f)\n',median([templates.template.r_like]),min([templates.template.r_like]),max([templates.template.r_like])))

            index=get_cell_indices(datarun,templates.cell_ids);
            templates.sptimes=[];  %all spike times of all 40 neurons
            templates.template_index=[];  %  the index of these 40 neurons as related to spike tme
            for i=1:length(templates.cell_ids)
                templates.sptimes=[templates.sptimes; datarun.spikes{index(i)}];
                templates.template_index=[templates.template_index; ones(size(datarun.spikes{index(i)}))*i]; 
            end
  
        end

       %all 40/100 whatever are selected... now wer are going to make
       %tmeplates that we will subtract
        
        if 1 %fast ei calculation - make templates use ei
            nlPoints=30;
            nrPoints=90;
            spike_nr_template=5000;
            sig_thr1=30;%pick elec for templates
            sig_thr2=0.5;%pick elec for fit
            sig_thr3=2;%mask in time
 
            shift_bin=.1;
            tshift=[-6:6];
            
            if ~exist([output_data_path '/temp/temp.ei'])  %if we don't already have the electrical images file
                display('calculate EIs')
                % save temp.neuron
                    old_neurons = edu.ucsc.neurobiology.vision.io.NeuronFile(datarun.names.rrs_neurons_path);
                    header = old_neurons.getHeader;
                    triggers = old_neurons.getTTLTimes;
                    old_neurons.close;  %i guess close it?

                    index=get_cell_indices(datarun,templates.cell_ids);
                    clear spikes
                    for i=1:length(index)
                        spikes{i}=int32(datarun.spikes{index(i)}*samplingrate); %convert from secs to time units of firing events
                    end

                    unix(['mkdir ' output_data_path '/temp']); %make directory
                    unix(['cp ' output_data_path '/*.globals ' output_data_path '/temp/temp.globals']);   % more direct unix commands on creating folders/files
                    save_neurons([output_data_path '/temp/temp.neurons'], header, triggers, spikes, datarun.cell_ids(index), datarun.channels(index));

                    %now carrying out the EI program in vision to find electrical image
                    command=sprintf('/snle/lab/Development/scripts/vision-auto-ei-only-64 %s %s/temp %d %d %d',raw_data_path,output_data_path, nlPoints, nrPoints, spike_nr_template);
                    unix(command);
            end
            
            %LOAD EI
            clear dataset
            dataset.names.rrs_ei_path=[output_data_path '/temp/temp.ei'];
            dataset.names.rrs_globals_path=[output_data_path '/temp/temp.globals'];
            dataset.names.rrs_neurons_path=[output_data_path '/temp/temp.neurons'];
            opt=struct('verbose',1,'load_neurons',1,'load_ei',1);
            dataset=load_data(dataset,opt);      %EI's are now calculated,  now we need to choose which part of the EI is real, and which part is garbage

            
            %We now have real EI's from which we make our deletion template
            T=text_waitbar('make templates');
            
            for j=1:length(templates.cell_ids)  %the top 40 neurons
                
                index=get_cell_indices(dataset,templates.cell_ids(j));   %now we have an actual single index
                ei=dataset.ei.eis{index};  %the actual  EI   (519 by 121   wehere 519 is the electrodes 121 is time 6 ms of time)
                ei(:,2:end)=ei(:,1:end-1); %hack match old code
                
                %thresholds
                [t1]=(max(abs(ei),[],2));
                [junk t2]=sort(t1);  %t2 is electrode index
                t3=ei(t2(8:end-1),1:20);  % first seven electrodes are blank recording only zeroes
                tstd=robust_std(t3(:));   %only using first msec to generate some sort of standard deviation  (maybe to characterize "noise" of EI)
                
                
  %              %sig_elec indexes the elctrodes from which we'll subtract
                sig_elec=find(t1>tstd*sig_thr1);     
                if length(sig_elec)<2
                    sig_elec=t2(end);  %then just use t
                end
                
  %              
                %i don't know what this does... will worry later
                [ttt1 ttt2]=sort(t1);
                if ttt1(end-1)>ttt1(end)*sig_thr2   %we say half for now
                    sig_elec_fit=ttt2(end-1:end);
                else
                    sig_elec_fit=ttt2(end);
                end 
                %so we either make not of the main spike... or the main
                %spike and the next largest neightbor if its more than half
            
                
                

                
                ttemplate=ei(sig_elec,:);  % on 519 scale raw subtraction... but need to ake into account aliasing the signal from sampling error
                
                %mask in time
                for i=1:length(sig_elec)
                    
                    %does not work reliable ttemplate(i,:)=ttemplate(i,:)-mean(ttemplate(i,1:15));
   %                 
                    t=circshift(conv(abs(ttemplate(i,:)),ones(1,20)/20,'same')',-10);   %  smoothed signal magnitude, smoothed over 1 msec
                    tt=find(t<tstd*sig_thr3);   %kill off small part of smoothed signal   key param... sig_thr3
                    ttemplate(i,tt)=0;
                    
                    if 0%plotting the individual electrode signal for a given neuron whose ei we are trying to erase
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
                templates.template(j).sig_elec_fit=sig_elec_fit+1; %shift for trigger that's in the data
                for i=1:length(templates.template(j).sig_elec_fit)   %second column of sig_elec_fit is the sig_elec index of sig_elec_fit electrodes
                    templates.template(j).sig_elec_fit(i,2)=find(templates.template(j).sig_elec==templates.template(j).sig_elec_fit(i,1));  %the index of the guy in sig elec
                end
                
                %shift templates
                x=[1:size(ttemplate,2)];  %number of columns... aka time steps   vector of it..   for now 1 2 3....121
                xi=[1:shift_bin:size(ttemplate,2)];    %vector  now with intervals one tenth the size
                ti=interp1(x,ttemplate',xi,'cubic');
                for i=1:length(tshift)
                    tti=circshift(ti,tshift(i));  %shifting in time, but by a tenth in time
                    ttii=interp1(xi,tti,x,'cubic');  %don't really need to interpolate... just picking off pts that have been shifted ever so slightly
                    templates.template(j).shift(i).mean=int16(ttii(1:end-1,:));     %is this the main info?   what's all this shifting do anyway?
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
                        plot(positions(templates.template(j).electrode,1),positions(templates.template(j).electrode,2),'y*')
                        plot(positions(sig_elec_fit,1),positions(sig_elec_fit,2),'g.')
                    subplot(2,1,2)
                        plot(templates.template(j).shift(1).mean)
                        
                    %pause
                end 

                if debug %plots
                    clf;
                    subplot(1,2,1)
                        positions = electrode_positions(array);
                        plot_ei_(ei,positions,0,'scale',5,'flipud',0,'cutoff',0.0,'sub_scale',.0) 
                        hold on
                        plot(positions(sig_elec,1),positions(sig_elec,2),'r.')
                        plot(positions(templates.template(j).electrode,1),positions(templates.template(j).electrode,2),'y*')
                        plot(positions(sig_elec_fit,1),positions(sig_elec_fit,2),'g.')
                    subplot(1,2,2)
                        plot(templates.template(j).shift(1).mean)
                    %saveas(gcf,sprintf('%s/template_%03d_%03d_%04d.jpg',output_data_path,j,templates.template(j).electrode,templates.template(j).cell_id),'jpg')   
                    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 20 10]); 
                    print('-djpeg',sprintf('%s/template_%03d_%03d_%04d.jpg',output_data_path,j,templates.template(j).electrode,templates.template(j).cell_id), '-r100');
                end
                
                T=text_waitbar(T,j/length(templates.cell_ids));
            end
            
            numsigelec=zeros(size(templates.cell_ids));
            for i=1:length(templates.cell_ids)
                numsigelec(i)=length(templates.template(i).sig_elec);  %the number of sig_elec it hits
            end 
            
            display(sprintf('median number of significant electrode per template: %d (min:%d max:%d)',median(numsigelec),min(numsigelec),max(numsigelec)))   

        end
        
 
        
 
        if 1 %save 
            display(['save ' output_data_path '/templates.mat']);
            save([output_data_path '/templates.mat'],'templates');
        end
        
        %test [s1 s2]=sort(templates.sptimes);templates.sptimes=templates.sptimes(s2);templates.template_index=templates.template_index(s2);
        
        
        
        
        if 1 %substract in raw data
            verbose=0;
            
            %re-use old header
            if ~exist('rdfh')   %raw data file header
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
                        %section where you see the actual subtraction
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
                        %pause
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
    data_path='/Volumes/Stream-phoenix/Data/Greschner/2011-08-04-4/data000';
    %data_path='/Volumes/Stream-phoenix/MacintoshHD/Data/Greschner/2011-08-04-4/data000';
    
    nr=99;
    
    import('edu.ucsc.neurobiology.vision.matlab.*');
    import('edu.ucsc.neurobiology.vision.io.*');
   
    tnames=dir([data_path '/serial_*']);
    
    sampling_frequency=edu.ucsc.neurobiology.vision.util.VisionParams.samplesPerMillisecond * 1000;
    
    new_neurons_file = [data_path '/serial.neurons'];
    
    %load data
    if 1
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
    if 1 
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
    
    if 0
        clf; set(gcf, 'color', 'white');

        ttlike=tlike;
        ttlike(find(~isfinite(ttlike)))=max(ttlike(find(isfinite(ttlike))));
        ttlike=ttlike/max(ttlike);

        subplot(2,1,1)
            hold on
            plot(trate/datarun.duration)
            plot(ttlike*max(trate/datarun.duration),'r');
            ylabel('rate')

        subplot(2,1,2)
            hold on
            plot(tamp)
            plot(ttlike*double(max(tamp)),'r');
            ylabel('amp')
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



            %{   
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
            %}
            %{
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
                
            %}  

                %{
                std of likelihood
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
           %}  




%{
% serial_01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

28-Sep-2011 14:20:08

input: /Volumes/Ground/Data/Akheitman/PlayingAround/data000
output: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_01/

use whitened covariance matrix
Creating output directory: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_01/

Running: Raw Data Noise Evaluation, at Wed Sep 28 14:20:09 PDT 2011
================= Progress Bar ===================
**************************************************
Raw Data Noise Evaluation done at Wed Sep 28 14:20:55 PDT 2011. Took 45.3sec.


Running: Spike Finding, at Wed Sep 28 14:20:55 PDT 2011
================= Progress Bar ===================
**************************************************
Spike Finding done at Wed Sep 28 15:04:52 PDT 2011. Took 2,636.9sec.


Running: Noise Whitened Covariances, at Wed Sep 28 15:04:52 PDT 2011
================= Progress Bar ===================
*************************************************
Electrode(s) [1, 130, 259, 260, 389, 390, 519] appear to be bad (had no spikes or singlular covariance matrices).

Noise Whitened Covariances done at Wed Sep 28 15:05:36 PDT 2011. Took 44.5sec.


Running: PCA Neuron Finding: Projections, at Wed Sep 28 15:05:36 PDT 2011
================= Progress Bar ===================
*************************************************
PCA Neuron Finding: Projections done at Wed Sep 28 15:12:49 PDT 2011. Took 432.7sec.


Running: PCA Neuron Finding: Clustering, at Wed Sep 28 15:12:49 PDT 2011
================= Progress Bar ===================
*************************************************
PCA Neuron Finding: Clustering done at Wed Sep 28 15:17:50 PDT 2011. Took 301.4sec.


Running: Neuron Cleaning, at Wed Sep 28 15:17:50 PDT 2011
Raw Neurons: 4076
After spike count cut: 4075
After contamination cut: 1833
Unique Neurons: 582
Neuron Cleaning done at Wed Sep 28 15:19:21 PDT 2011. Took 90.8sec.

/Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_01//serial_01.globals

Running: Make Parameters File, at Wed Sep 28 15:19:21 PDT 2011
================= Progress Bar ===================
*************************************Auto: 63 %
Contamination Index: 36 %

Make Parameters File done at Wed Sep 28 15:19:22 PDT 2011. Took .7sec.


The whole analysis took: 59.2 min.
load neurons file: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_01//serial_01.neurons-raw
Examining 4077 cells (RRS v.32) ... extracted 4076 cells.

|-find clusters----------------------------------|
..................................................

select clusters...
cluster found: 1472
cluster selected: 40
median cluster likelihood: 159.70 (min:109.15 max:475.02)

calculate EIs
Saving 40 neurons ... done.
Parameter 'Dataset Folder' is set to: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_01//temp
Parameter 'Raw Data File' is set to: /Volumes/Ground/Data/Akheitman/PlayingAround/data000
Parameter 'Mean Time Constant' is set to: 0.01
Parameter 'Left Samples' is set to: 30
Parameter 'Right Samples' is set to: 90
Parameter 'Spikes To Average' is set to: 5000
Parameter 'nThreads' is set to: 4
Running: Electrophysiological Imaging Fast

Running: Electrophysiological Imaging Fast, at Wed Sep 28 15:22:20 PDT 2011
================= Progress Bar ===================
**************************************************
Electrophysiological Imaging Fast done at Wed Sep 28 15:25:49 PDT 2011. Took 209.1sec.

load neurons file: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_01//temp/temp.neurons
Examining 41 cells (RRS v.32) ... extracted 40 cells.
load ei file: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_01//temp/temp.ei
!!!!!! ARRAY ID WAS LOADED FROM GLOBALS: 1601 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!   WARNING: array IDs have not been set yet for values >= 1500   !!!!!!!!!!!
!!!!!!!!   In addition, this is an unknown > 1500 array id;              !!!!!!!!!!!
!!!!!!!!       transforms may be broken                                  !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

|-make templates---------------------------------|
..................................................

median significant electrodes: 3.250000e+01 (min:6 max:122)
save /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_01//templates.mat

---> Saving Raw Data in 1 file(s):
Location 1: //Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_01//bin//2011-08-04-4/data000.bin

|-process bin------------------------------------|
..................................................

rm: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_-1/bin/: No such file or directory
Elapsed time is 5040.087383 seconds.


% serial_02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

28-Sep-2011 15:44:09

input: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_01//bin/data999/
output: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_02/

use whitened covariance matrix
Creating output directory: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_02/

Running: Raw Data Noise Evaluation, at Wed Sep 28 15:44:40 PDT 2011
================= Progress Bar ===================
**************************************************
Raw Data Noise Evaluation done at Wed Sep 28 15:45:24 PDT 2011. Took 43.6sec.


Running: Spike Finding, at Wed Sep 28 15:45:24 PDT 2011
================= Progress Bar ===================
**************************************************
Spike Finding done at Wed Sep 28 16:28:56 PDT 2011. Took 2,611.8sec.


Running: Noise Whitened Covariances, at Wed Sep 28 16:28:56 PDT 2011
================= Progress Bar ===================
*************************************************
Electrode(s) [1, 130, 259, 260, 389, 390, 519] appear to be bad (had no spikes or singlular covariance matrices).

Noise Whitened Covariances done at Wed Sep 28 16:29:44 PDT 2011. Took 47.4sec.


Running: PCA Neuron Finding: Projections, at Wed Sep 28 16:29:44 PDT 2011
================= Progress Bar ===================
*************************************************
Warning: electrode 245 has zero spikes, check thresholds, they may be too high
PCA Neuron Finding: Projections done at Wed Sep 28 16:36:41 PDT 2011. Took 417.6sec.


Running: PCA Neuron Finding: Clustering, at Wed Sep 28 16:36:41 PDT 2011
================= Progress Bar ===================
*************************************************
PCA Neuron Finding: Clustering done at Wed Sep 28 16:42:08 PDT 2011. Took 326.5sec.


Running: Neuron Cleaning, at Wed Sep 28 16:42:08 PDT 2011
Raw Neurons: 4067
After spike count cut: 4066
After contamination cut: 1822
Unique Neurons: 584
Neuron Cleaning done at Wed Sep 28 16:43:36 PDT 2011. Took 88.3sec.

/Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_02//serial_02.globals

Running: Make Parameters File, at Wed Sep 28 16:43:36 PDT 2011
================= Progress Bar ===================
*************************************Auto: 62 %
Contamination Index: 37 %

Make Parameters File done at Wed Sep 28 16:43:37 PDT 2011. Took .7sec.


The whole analysis took: 59.0 min.
load neurons file: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_02//serial_02.neurons-raw
Examining 4068 cells (RRS v.32) ... extracted 4067 cells.

|-find clusters----------------------------------|
..................................................

select clusters...
cluster found: 1453
cluster selected: 40
median cluster likelihood: 127.61 (min:85.29 max:Inf)

calculate EIs
Saving 40 neurons ... done.
Parameter 'Dataset Folder' is set to: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_02//temp
Parameter 'Raw Data File' is set to: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_01//bin/data999/
Parameter 'Mean Time Constant' is set to: 0.01
Parameter 'Left Samples' is set to: 30
Parameter 'Right Samples' is set to: 90
Parameter 'Spikes To Average' is set to: 5000
Parameter 'nThreads' is set to: 4
Running: Electrophysiological Imaging Fast

Running: Electrophysiological Imaging Fast, at Wed Sep 28 16:47:15 PDT 2011
================= Progress Bar ===================
**************************************************
Electrophysiological Imaging Fast done at Wed Sep 28 16:51:26 PDT 2011. Took 250.7sec.

load neurons file: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_02//temp/temp.neurons
Examining 41 cells (RRS v.32) ... extracted 40 cells.
load ei file: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_02//temp/temp.ei
!!!!!! ARRAY ID WAS LOADED FROM GLOBALS: 1601 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!   WARNING: array IDs have not been set yet for values >= 1500   !!!!!!!!!!!
!!!!!!!!   In addition, this is an unknown > 1500 array id;              !!!!!!!!!!!
!!!!!!!!       transforms may be broken                                  !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

|-make templates---------------------------------|
..................................................

median significant electrodes: 38 (min:9 max:76)
save /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_02//templates.mat

---> Saving Raw Data in 1 file(s):
Location 1: //Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_02//bin//2011-08-04-4/data000.bin

|-process bin------------------------------------|
..................................................

rm: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_00/bin/: No such file or directory
Elapsed time is 5605.935681 seconds.


% serial_03 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

28-Sep-2011 17:17:41

input: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_02//bin/data999/
output: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_03/

use whitened covariance matrix
Creating output directory: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_03/

Running: Raw Data Noise Evaluation, at Wed Sep 28 17:18:03 PDT 2011
================= Progress Bar ===================
**************************************************
Raw Data Noise Evaluation done at Wed Sep 28 17:18:45 PDT 2011. Took 41.9sec.


Running: Spike Finding, at Wed Sep 28 17:18:45 PDT 2011
================= Progress Bar ===================
**************************************************
Spike Finding done at Wed Sep 28 18:01:39 PDT 2011. Took 2,574.5sec.


Running: Noise Whitened Covariances, at Wed Sep 28 18:01:39 PDT 2011
================= Progress Bar ===================
*************************************************
Electrode(s) [1, 130, 259, 260, 389, 390, 519] appear to be bad (had no spikes or singlular covariance matrices).

Noise Whitened Covariances done at Wed Sep 28 18:02:24 PDT 2011. Took 44.6sec.


Running: PCA Neuron Finding: Projections, at Wed Sep 28 18:02:24 PDT 2011
================= Progress Bar ===================
*************************************************
Warning: electrode 245 has zero spikes, check thresholds, they may be too high
PCA Neuron Finding: Projections done at Wed Sep 28 18:09:54 PDT 2011. Took 449.6sec.


Running: PCA Neuron Finding: Clustering, at Wed Sep 28 18:09:54 PDT 2011
================= Progress Bar ===================
*************************************************
PCA Neuron Finding: Clustering done at Wed Sep 28 18:14:24 PDT 2011. Took 270.6sec.


Running: Neuron Cleaning, at Wed Sep 28 18:14:24 PDT 2011
Raw Neurons: 4071
After spike count cut: 4071
After contamination cut: 1840
Unique Neurons: 571
Neuron Cleaning done at Wed Sep 28 18:15:52 PDT 2011. Took 88.3sec.

/Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_03//serial_03.globals

Running: Make Parameters File, at Wed Sep 28 18:15:52 PDT 2011
================= Progress Bar ===================
************************************Auto: 61 %
Contamination Index: 38 %

Make Parameters File done at Wed Sep 28 18:15:53 PDT 2011. Took .7sec.


The whole analysis took: 57.9 min.
load neurons file: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_03//serial_03.neurons-raw
Examining 4072 cells (RRS v.32) ... extracted 4071 cells.

|-find clusters----------------------------------|
..................................................

select clusters...
cluster found: 1447
cluster selected: 40
median cluster likelihood: 94.73 (min:68.11 max:299.82)

calculate EIs
Saving 40 neurons ... done.
Parameter 'Dataset Folder' is set to: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_03//temp
Parameter 'Raw Data File' is set to: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_02//bin/data999/
Parameter 'Mean Time Constant' is set to: 0.01
Parameter 'Left Samples' is set to: 30
Parameter 'Right Samples' is set to: 90
Parameter 'Spikes To Average' is set to: 5000
Parameter 'nThreads' is set to: 4
Running: Electrophysiological Imaging Fast

Running: Electrophysiological Imaging Fast, at Wed Sep 28 18:19:31 PDT 2011
================= Progress Bar ===================
**************************************************
Electrophysiological Imaging Fast done at Wed Sep 28 18:23:51 PDT 2011. Took 259.8sec.

load neurons file: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_03//temp/temp.neurons
Examining 41 cells (RRS v.32) ... extracted 40 cells.
load ei file: /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_03//temp/temp.ei
!!!!!! ARRAY ID WAS LOADED FROM GLOBALS: 1601 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!   WARNING: array IDs have not been set yet for values >= 1500   !!!!!!!!!!!
!!!!!!!!   In addition, this is an unknown > 1500 array id;              !!!!!!!!!!!
!!!!!!!!       transforms may be broken                                  !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

|-make templates---------------------------------|
..................................................

median significant electrodes: 3.150000e+01 (min:9 max:79)
save /Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_03//templates.mat

---> Saving Raw Data in 1 file(s):
Location 1: //Volumes/Ground/Data/Akheitman/PlayingAround/data000/serial_03//bin//2011-08-04-4/data000.bin

|-process bin------------------------------------|
..................................................

Elapsed time is 5336.319896 seconds.
%}










