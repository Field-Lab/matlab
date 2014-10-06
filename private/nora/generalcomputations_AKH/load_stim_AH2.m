function datarun = load_stim_AH2(datarun,varargin) 
%load obvious stimulus files 
%greschner

% AH comments for myself
% .params is kind of a meaningless entry  just the list of start frames
% .combinations is a list of the first group of start frames
% .trials is the full list of start frames
% .trails_list is the full list of start frames when translated to the index
% of the combinations



% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('ignor_field', []);%eg {'SEED'}
    p.addParamValue('change_entry', []);% struct('field','SEED','find',{11112,11113},'replace',{11111,11111})
    p.addParamValue('find_trigger', 0);
    p.addParamValue('trigger_iti_thr', 2);
    p.addParamValue('plot', 0);
    p.addParamValue('correction_incomplet_run', 1);

    p.parse(varargin{:});
    params = p.Results;
    
  

clear s
            fid=fopen(datarun.names.stimulus_path);
            
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
                        
                        if ~isempty(str2num((l{i+1}))) | isequal(l{i+1},'0')
                            s.trials(jj).(t)=str2num(l{i+1});
                        else
                            s.trials(jj).(t)=(l{i+1});
                        end
                        n=2;
                        while i+n<=length(l) & isempty(strfind(l{i+n},':'))
                            if ~isempty(str2num((l{i+n})))
                                s.trials(jj).(t)=[s.trials(jj).(t) str2num(l{i+n})];
                            else
                                s.trials(jj).(t)=[s.trials(jj).(t) [ ] (l{i+n})];
                            end
                            n=n+1;
                        end
                    end
            end
            
            
  
           
            %combinations
           
            %ignor
            if ~isempty(params.ignor_field)
                save_trails=s.trials;
                if iscell(params.ignor_field)
                    for jj=1:length(params.ignor_field)
                        for j=1:length(s.trials)
                            s.trials(j).(params.ignor_field{jj})=-i;
                        end
                    end
                else
                    for j=1:length(s.trials)
                        s.trials(j).(params.ignor_field)=-i;
                    end    
                end
            end
                
            %change_entry   
            if ~isempty(params.change_entry)
                for j=1:length(s.trials)
                    for jj=1:length(params.change_entry)
                        if ~isequal(s.trials(j).(params.change_entry(jj).field),params.change_entry(jj).find)
                            s.trials(j).(params.change_entry(jj).field)=params.change_entry(jj).replace;
                        end
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
   
            %ignor
            if ~isempty(params.ignor_field)
                s.trials=save_trails;
            end
            
            
  
if mod(length(s.trials),length(s.combinations))
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('load_stim: uneven number of repetions');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
end              
                       
disp(sprintf('load_stim: %d stimulus combinations   %d trials   %d repetitions',length(s.combinations),length(s.trials),s.repetitions));
disp(s.params)





%trigger%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i have no idea what is going here  AH %
if params.find_trigger
    %finds change in inter trigger interval
        tt=datarun.triggers(2:end)-datarun.triggers(1:end-1);
        thr=std(tt)*params.trigger_iti_thr;
        ttt=datarun.triggers(1);
        alt=tt(1);
        for ii=2:length(tt)
            if abs(tt(ii)-alt)>thr
                ttt=[ttt datarun.triggers(ii+1)];
                if ii<length(tt)
                    ii=ii+1;
                end
            end
            alt=tt(ii);
        end
        triggers=ttt;


        if params.plot
            plot(datarun.triggers,'.');
            hold on
            plot([ones(size(triggers)); ones(size(triggers))*length(datarun.triggers)],[triggers; triggers],'r');
        end
else
    triggers=datarun.triggers;
end

%missing trigger
if length(triggers)~=length(s.trials)
            disp(sprintf('!!!!!!!!!!!!\n load_stim: ERROR %d triggers - %d stimuli\n!!!!!!!!!!!!',length(triggers), length(s.trials)));

            if params.correction_incomplet_run
                trepeats=floor(length(triggers)/length(s.combinations));
                disp(sprintf('\n asume run ended early  delete last repetions: %d rep -> %d rep\n',s.repetitions, trepeats)); 
                s.repetitions=trepeats;
                s.trials=s.trials(1:s.repetitions*length(s.combinations));
                s.trial_list=s.trial_list(1:s.repetitions*length(s.combinations));
                triggers=triggers(1:s.repetitions*length(s.combinations));

                if params.plot
                    plot([ones(size(triggers)); ones(size(triggers))*length(datarun.triggers)],[triggers; triggers],'g');
                    hold off
                end  

            end
end
s.triggers=triggers;
     
     
     
 t=fieldnames(s);
 for i=1:length(t)
    datarun.stimulus.(t{i})=getfield(s,t{i});
 end
               

 







    







