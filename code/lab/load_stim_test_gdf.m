function datarun = load_stim(datarun,varargin) 
%load obvious stimulus files 
%greschner
%sravi - 03-2014 - modified function to be compatible with newer version of
                   %matlab (textscan function not reading the string in stimulus file - use fileread)

% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('ignor_field', []);%eg {'SEED'}
    p.addParamValue('change_entry', []);% struct('field','SEED','find',{11112,11113},'replace',{11111,11111})
    p.addParamValue('find_trigger', 1);
    p.addParamValue('trigger_iti_thr', 2);
    p.addParamValue('plot', 0);
    p.addParamValue('correction_incomplet_run', 1);
    p.addParamValue('user_defined_trigger_interval', [], @isnumeric); % Caution, this might not work if frames were dropped
    p.addParamValue('user_defined_trigger_set', [], @isnumeric);
    
    p.parse(varargin{:});
    params = p.Results;
    
  

clear s

            %textscan modification here 
            text = fileread(datarun.names.stimulus_path); %returns contents of file as matlab string
            if(isspace(text(end)))
                text(end) = [];  %removes space character at the end of stimulus files
            end
            file=textscan(text, '%s'); %read the string
            l=file{1};

            %parse first line, a.k.a header parameters
            t=1;
            j=2;
            while t
                if ~isempty(strfind(l{j},')'));
                        t=t-1;
                end
                if ~isempty(strfind(l{j},'))'));
                        t=t-1;
                end
                if ~isempty(strfind(l{j},'('));
                        t=t+1;
                end
                j=j+1;
            end
            
            
            for i=1:length(l)
                l{i}=strrep(l{i},'(','');
                l{i}=strrep(l{i},')','');
                l{i}=strrep(l{i},'#','');
            end
           
            jj=0;
            for i=j:length(l)
                    if ~isempty(strfind(l{i},':'));
                        t=strrep(l{i},':','');
                        t=strrep(t,'-','_');
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
                    %gdf modified this 'else' statement: Martin's old code
                    %is commented out below. 2013-04-22
                    % get parameter values
                    temp_vals = zeros(1,length(s.combinations));
                    for cnd = 1:length(s.combinations);
                        temp_vals(cnd) = s.combinations(cnd).(t{j})(1);
                    end
                    temp_vals = sort(temp_vals, 'ascend');
                    s.params.(t{j}) = temp_vals;
                        
                    % martin's old code -- gdf doesn't think it does the
                    % right thing
%                     s.params.(t{j})=s.combinations(1).(t{j});
%                     tt=0;
%                     for jj=2:length(s.combinations)
%                         for jjj=1:size(s.params.(t{j}),1)  
%                             if all(s.combinations(jj).(t{j})(jjj,:)==s.params.(t{j}))
%                                 tt=1;
%                             end
%                         end
%                         if ~tt
%                             s.combinations(jj).(t{j})(jj,:)=s.params.(t{j});
%                         end
%                     end
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

if ~isempty(p.Results.user_defined_trigger_set);
    triggers = datarun.triggers(p.Results.user_defined_trigger_set);
elseif ~isempty(p.Results.user_defined_trigger_interval)
    % define new stimulus trigger interval
    stim_interval = p.Results.user_defined_trigger_interval; % units are seconds
    % find the triggers that are closest to occuring subsequent to these intervals
    stim_trigger_indices = find(mod(datarun.triggers,stim_interval) < 0.1);
    % set the triggers
    triggers = datarun.triggers(stim_trigger_indices)';

    %trigger%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif params.find_trigger && isempty(p.Results.user_defined_trigger_interval)
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
               
