function ds=DS_load_stim_pair_consecutive_repeats(file, data) 

%file=stim_file;
%clear s stimulus ds

%load stimulus
    
            triggers=data.triggers;
            %tempperiod=data.ds.tempperiod;
     
            fid=fopen(file);
            
            file=textscan(fid, '%s');
            fclose(fid);
            l=file{1};

            for i=1:length(l)
                l{i}=strrep(l{i},'(','');
                l{i}=strrep(l{i},')','');
                l{i}=strrep(l{i},':','');
                l{i}=strrep(l{i},'#','');
            end
    
            
            i=1;
            while isempty(strfind(l{i},'SPATIAL-PERIOD'))
               i=i+1; 
            end      
            n=0;
            for i=i:length(l)
                if ~isempty(strfind(l{i},'SPATIAL-PERIOD'));
                    n=n+1;
                    s.spatial(n)=str2num(l{i+1});  
                end
                if ~isempty(strfind(l{i},'TEMPORAL-PERIOD'));
                    s.temporal(n)=str2num(l{i+1});
                end
                if ~isempty(strfind(l{i},'DIRECTION'));
                    s.direction(n)=str2num(l{i+1});
                end
                if ~isempty(strfind(l{i},'RGB'));
                    s.rgb(n)=str2num(l{i+1});
                end
            end
            
            
               %bug fix stim 060705-0/007
               if 0
                    sh=ones(1,8);
                    for i=1:length(s.spatial)
                        if isequal(s.spatial(i),256) & isequal(s.temporal(i),64)
                            if sh(s.direction(i)/45+1)
                                sh(s.direction(i)/45+1)=0;
                            else
                                %s.temporal(i)=1256;
                                s.spatial(i)=1256;
                                sh(s.direction(i)/45+1)=1;
                            end  
                        end
                    end 
               end  
               
               
            %how many
            t=[];tt=[];ttt=[];tttt=[];
            for i=1:length(s.spatial)
                if ~ismember(s.spatial(i), t)
                    t=[t s.spatial(i)];
                end
                if ~ismember(s.temporal(i), tt)
                    tt=[tt s.temporal(i)];
                end
                if ~ismember(s.direction(i), ttt)
                    ttt=[ttt s.direction(i)];
                end
                if isfield(s, 'rgb')
                    if ~ismember(s.rgb(i), tttt)
                        tttt=[tttt s.rgb(i)];
                    end
                end
            end
            spatial=sort(t);
            temporal=sort(tt);
            direction=sort(ttt);
            if isfield(s, 'rgb')
                rgb=sort(tttt);
            else
                for i=1:length(l)
                    if ~isempty(strfind(l{i},'RGB'))
                        rgb=str2num(l{i+1});
                    end
                end
                for i=1:length(s.direction)
                    s.rgb(i)=rgb;
                end
            end
            
            %how may pairs  assumption if color, then in every pair 
            n=1;
            for i=1:length(spatial)
                t1=find(s.spatial==spatial(i));
                t2=[];
                for ii=1:length(t1)
                    if ~ismember(s.temporal(t1(ii)), t2)
                        t2=[t2 s.temporal(t1(ii))];
                        for iii=1:length(rgb)
                            ds.stimulus(n).spatial=s.spatial(t1(ii));
                            ds.stimulus(n).temporal=s.temporal(t1(ii));
                            ds.stimulus(n).rgb=rgb(iii);
                            n=n+1;
                        end
                    end
                end  
            end
            

                            
            Directions=length(direction);
            Repeats=length(s.direction)/length(ds.stimulus)/Directions;
            
            ds.repeats=Repeats;
            ds.number_directions=Directions;
            ds.directions=direction;
            
            
            tt=triggers(2:end)-triggers(1:end-1);
            ttt=triggers(1);
            alt=tt(1);
            for ii=2:length(tt)
                if (tt(ii)-alt)^2>.001
                    ttt=[ttt triggers(ii+1)];
                    if ii<length(tt)
                        ii=ii+1;
                    end
                end
                alt=tt(ii);
            end
            if length(ttt)~=ds.repeats*ds.number_directions*length(ds.stimulus)
                disp(sprintf('!!!!!!!!!!!!\n!!\n!!  DS_load_stim_pair: ERROR\n!!\n!!!!!!!!!!!!'));
            end

            
%obvious to matlab: screen 0 top 
            %s.direction=-s.direction+90;
            %s.direction(find(s.direction<0))=360+s.direction(find(s.direction<0));
            
           
            for j=1:length(ds.stimulus)
                for jj=1:length(direction)
                    s1=find( s.spatial==ds.stimulus(j).spatial & s.temporal==ds.stimulus(j).temporal & s.rgb==ds.stimulus(j).rgb & s.direction==direction(jj) );
                    ds.stimulus(j).tr(jj,:)=ttt(sort(s1));
                end
            end
            



            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
