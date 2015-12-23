 if 0 %load
    if 1
        clear datarun
        datarun{1}.names.rrs_params_path='/snle/lab/Experiments/Array/Analysis/2011-09-14-1/data000/data000.params';
        datarun{2}.names.rrs_neurons_path='/snle/lab/Experiments/Array/Analysis/2011-09-14-1/data002-from-data000/data002-from-data000.neurons';
        datarun{2}.names.stimulus_path='/snle/lab/Experiments/Array/Analysis/2011-09-14-1/data002-from-data000/s02';
    end

    opt=struct('verbose',1,'load_params',1,'load_neurons',1);
    datarun=load_data(datarun,opt);
    if length(datarun)>1
        datarun=map_cell_types(datarun, struct('verbose',true));  
        datarun{2}.ds=DS_load_stim_pair_consecutive_repeats(datarun{2}.names.stimulus_path, datarun{2}); 
    else
        datarun{1}.ds=DS_load_stim_pair_consecutive_repeats(datarun{1}.names.stimulus_path, datarun{1});
        %datarun{1}=DS_load_stim_pair(datarun{1});
    end
end    


if 1 %plot table both screen sort
   nr=2; 
   clf; set(gcf, 'color', 'white')
   
   sp=2;re=1;
   
   start=0.5; stop=8;  
   index=[1:length(datarun{nr}.cell_ids)];
   
   if 1 %compute 
       clear list
       for j=1:length(datarun{nr}.ds.stimulus)
           for i=1:length(index)
                t1=zeros(1,datarun{nr}.ds.number_directions);
                for ii=1:datarun{nr}.ds.number_directions
                    for iii=1:datarun{nr}.ds.repeats 
                        h=datarun{nr}.spikes{index(i)}'-datarun{nr}.ds.stimulus(1).tr(ii,iii);  
                        t1(ii)=t1(ii)+length(find(h>=start & h<=stop));
                    end
                end
                list(i).stimulus(j).average=t1;
                list(i).stimulus(j).x=sum(cos(datarun{nr}.ds.directions/180*pi).*t1)/sum(t1);
                list(i).stimulus(j).y=sum(sin(datarun{nr}.ds.directions/180*pi).*t1)/sum(t1);
                [list(i).stimulus(j).theta list(i).stimulus(j).rho]=cart2pol(list(i).stimulus(j).x,list(i).stimulus(j).y);     
           end
       end
       [junk,t]=sort(rho);
       index=index(flipud(t));
   end
   

   k=1; kmin=1; kmax=length(index); ha=loop_slider(k,kmin,kmax);
   while k
        k=round(get(ha,'Value')); 

        uu=index(k);

        for i=1:sp
            [junk,s]=sort([datarun{nr}.ds.stimulus((i-1)*re+1:i*re).temporal]);
            for ii=1:re

                subplot(sp*2,re,(i-1)*re*2+ii);
                ttr=reshape(datarun{nr}.ds.stimulus((i-1)*re+s(ii)).tr',1,datarun{nr}.ds.number_directions*datarun{nr}.ds.repeats);
                %if ii==1
                %    psth_raster(0,Length/4,datarun{nr}.spikes{uu}',ttr);
                %end
                if ii<=3
                    psth_raster(start,stop,datarun{nr}.spikes{uu}',ttr);
                end
                if ii>3
                    psth_raster(start,stop,datarun{nr}.spikes{uu}',ttr);
                end
                title(sprintf('spatial=%d  temporal=%d ',datarun{nr}.ds.stimulus((i-1)*re+s(ii)).spatial,datarun{nr}.ds.stimulus((i-1)*re+s(ii)).temporal));     

                subplot(sp*2,re,(i-1)*re*2+re+ii);
                clear tt
                for iii=1:datarun{nr}.ds.number_directions
                    ttt=0;
                    for iiii=1:datarun{nr}.ds.repeats 
                        h=datarun{nr}.spikes{uu}'-datarun{nr}.ds.stimulus((i-1)*re+s(ii)).tr(iii,iiii);  
                        ttt=ttt+length(find(h>=start & h<=stop));
                    end
                    tt(iii)=ttt;
                end
                hold off;
                mpolar((tt));

            end
        end
        title(sprintf('cell-id=%d',datarun{nr}.cell_ids(uu)));     

        uiwait;
	end
       
end


if 0
    old_file_name=datarun{1}.names.rrs_params_path;
    new_file_name=[old_file_name(1:end-7) '-DS.params'];
    col_names={'directions' 'averageGrating' 'xDS1' 'yDS1' 'magDS1' 'angDS1' 'xDS2' 'yDS2' 'magDS2' 'angDS2'};
    col_types={'DoubleArray' 'DoubleArray' 'Double' 'Double' 'Double' 'Double' 'DoubleArray' 'DoubleArray' 'Double' 'Double' 'Double' 'Double'};
    
    clear data
    for i=1:length(index)
        data{1,i}=datarun{nr}.ds.directions;
        
    end
    
    add2params(old_file_name,new_file_name,col_names,col_types,data)
    
end





























