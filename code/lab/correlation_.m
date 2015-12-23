function synchrony_index=correlation_(spikes_1, spikes_2, varargin) 
%returns matrix of cell_specification_1 x cell_specification_2 with the
%correlation strength log2(p_ab/(p_a*p_b))
%
% synchrony_index=correlation(datarun, cell_specification_1, [cell_specification_2], [params]) 
%
% defaults.bin = 10; 


% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('corrcoef', 1);%
    p.addParamValue('bin', 10);% 
    p.addParamValue('shift', 0);%if 2. spiketrain is shifted pos, shift is pos 
   
    % parse inputs
    p.parse(varargin{:});
    params = p.Results;
  
    params.shift=-params.shift/1000;
    
if exist('spikes_2','var') && ~isempty(spikes_2);
    duration=ceil(max([cell2mat(spikes_1); cell2mat(spikes_2)])); 
else
    duration=ceil(max([cell2mat(spikes_1)]));
end


%%%%jons paper style%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%not changed!!!!!
if ~params.corrcoef

    if isequal(index_1, index_2) %split for speed 

        for i=1:length(index_1)
            spiketrain_1=histc(datarun.spikes{index_1(i)},[0:params.bin/1000:datarun.duration-params.bin/1000]); 
            p_a=length(find(spiketrain_1))/length(spiketrain_1);
            for ii=i:length(index_2)

                spiketrain_2=histc(datarun.spikes{index_2(ii)},[0:params.bin/1000:datarun.duration-params.bin/1000]);
                p_ab=length(find(spiketrain_1 & spiketrain_2))/length(spiketrain_1);
                if p_ab
                    p_b=length(find(spiketrain_2))/length(spiketrain_2);
                    synchrony_index(index_1(i),index_2(ii))=log2(p_ab/(p_a*p_b));
                else
                    synchrony_index(index_1(i),index_2(ii))=0;
                end 
                synchrony_index(index_2(ii),index_1(i))=synchrony_index(index_1(i),index_2(ii));
                
            end       
        end

    else

        for i=1:length(index_1)
            spiketrain_1=histc(datarun.spikes{index_1(i)},[0:params.bin/1000:datarun.duration-params.bin/1000]); 
            p_a=length(find(spiketrain_1))/length(spiketrain_1);
            for ii=1:length(index_2)

                spiketrain_2=histc(datarun.spikes{index_2(ii)},[0:params.bin/1000:datarun.duration-params.bin/1000]);
                
                p_ab=length(find(spiketrain_1 & spiketrain_2))/length(spiketrain_1);
                if p_ab
                    p_b=length(find(spiketrain_2))/length(spiketrain_2);
                    synchrony_index(index_1(i),index_2(ii))=log2(p_ab/(p_a*p_b));
                else
                    synchrony_index(index_1(i),index_2(ii))=0;
                end

            end       
        end

    end 
end

    
%%%%corrcoef%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if params.corrcoef
     
    %split for speed vs memory 
    if ~ (exist('spikes_2','var') && ~isempty(spikes_2));
        try
            spiketrain=zeros(floor(duration*1000/params.bin),length(spikes_1));
            for i=1:length(spikes_1)
                spiketrain(:,i)=histc(spikes_1{i},[0:params.bin/1000:duration-params.bin/1000]);    
            end
            synchrony_index=corrcoef(spiketrain);
            fast=true;
        catch
            fast=false;
            spikes_2=spikes_1;
        end
    else
       try 
            spiketrain=zeros(floor(duration*1000/params.bin),length(spikes_1)+length(spikes_2));
            for i=1:length(spikes_1)
                spiketrain(:,i)=histc(spikes_1{i},[0:params.bin/1000:duration-params.bin/1000]);    
            end
            for i=1:length(spikes_2)
                spiketrain(:,i+length(spikes_1))=histc(spikes_2{i}+params.shift,[0:params.bin/1000:duration-params.bin/1000]);    
            end
            temp=corrcoef(spiketrain);
            synchrony_index=temp(1:length(spikes_1),length(spikes_1)+1:end);
            fast=true;
       catch
            fast=false;
       end
    end
    
    if ~fast
            T=text_waitbar('correlation');
            spiketrain=zeros(duration*1000/params.bin,2);
            synchrony_index=zeros(length(spikes_1),length(spikes_2));
            for i=1:length(spikes_1)
                spiketrain(:,1)=histc(spikes_1{i},[0:params.bin/1000:duration-params.bin/1000]); 
                for ii=1:length(spikes_2) 
                    spiketrain(:,2)=histc(spikes_2{ii},[0:params.bin/1000:duration-params.bin/1000]); 
                    tem=corrcoef(spiketrain);
                    temp(i,ii)=tem(1,2);
                end
                T=text_waitbar(T,i/length(spikes_1));
            end
            synchrony_index=temp;
    end
       
end
   
 
    
  
