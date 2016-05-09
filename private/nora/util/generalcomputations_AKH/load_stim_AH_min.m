% NOTE doesn't work on the altest matlab versions
% does work on Matlab 7.10 which is on Smokestack! 
function datarun = load_stim_AH_min(datarun)



            fid=fopen(datarun.names.stimulus_path);
            
          %  fid=fopen(dataruna.names.stimulus_path);
            
            file=textscan(fid, '%s');
            fclose(fid);
            l=file{1};

            %torse first line 
            t=1;
            j=2;
            % looping through each line .. liking for a '('   as in (:Start) 
            while t
                if ~isempty(strfind(l{j},')'));
                        t=t-1;
                end
                if ~isempty(strfind(l{j},'('));
                        t=t+1;
                end
                j=j+1;  
            end
            % j is a place holder for the start of the frame numbers 
            
            % just process l to get rid of all '(', ')', '#'
            for i=1:length(l)
                l{i}=strrep(l{i},'(','');
                l{i}=strrep(l{i},')','');
                l{i}=strrep(l{i},'#','');
            end
           
            jj=0; % jj is counter for the actual start frames
            % inroduce s of which
            %abuse of letter t
            % s.trials(start_frame)  seems to denote the start frame
            clear t
            for i=j:length(l)
                    if ~isempty(strfind(l{i},':'));
                        t=strrep(l{i},':',''); % get rid of semii-colon
                        t=strrep(t,'-','_'); % change - to _
                        if isequal(l{j},l{i})
                            jj=jj+1;
                        end
                        
                        if ~isempty(str2num((l{i+1}))) || isequal(l{i+1},'0')
                            s.trials(jj).(t)=str2num(l{i+1});
                        else
                            s.trials(jj).(t)=(l{i+1});
                        end
                        n=2;
                        while i+n<=length(l) && isempty(strfind(l{i+n},':'))
                            if ~isempty(str2num((l{i+n})))
                                s.trials(jj).(t)=[s.trials(jj).(t) str2num(l{i+n})];
                            else
                                s.trials(jj).(t)=[s.trials(jj).(t) [ ] (l{i+n})];
                            end
                            n=n+1;
                        end
                    end
            end
            
            
            
                        s.combinations(1)=s.trials(1);
            s.trial_list(1)=1;
            for j=2:length(s.trials)
                t=0;
                for jj=1:length(s.combinations)
                    if isequal(s.combinations(jj),s.trials(j))
                        s.trial_list(j)=jj;
                        t=jj;
                    end    
                end
                if ~t
                    s.combinations(length(s.combinations)+1)=s.trials(j);
                    s.trial_list(j)=length(s.combinations);
                end
            end

                       
            t=fieldnames(s.combinations);
            for j=1:length(t)
                if length(s.combinations(1).(t{j}))==1
                    s.params.(t{j})=sort(unique([s.combinations.(t{j})]));
                else
                    s.params.(t{j})=s.combinations(1).(t{j});
                    tt=0;
                    for jj=2:length(s.combinations)
                        for jjj=1:size(s.params.(t{j}),1)  
                            if all(s.combinations(jj).(t{j})(jjj,:)==s.params.(t{j}))
                                tt=1;
                            end
                        end
                        if ~tt
                            s.combinations(jj).(t{j})(jj,:)=s.params.(t{j});
                        end
                    end
                end
            end
            
            s.repetitions=length(s.trials)/length(s.combinations);
            
            
            
            
            
            
            
            triggers=datarun.triggers;
            
            
            
            
            
            
            s.triggers=triggers;
     
     
     
 t=fieldnames(s);
 for i=1:length(t)
    datarun.stimulus.(t{i})=getfield(s,t{i});
 end
 
 
 
end