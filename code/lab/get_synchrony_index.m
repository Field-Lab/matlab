function datarun=get_synchrony_index(datarun, cell_specification, varargin) 
%adds synchrony_index field:
%matrix of cell_index x cell_index with correlation coeff or correlation index log2(p_ab/(p_a*p_b))
% 
%
% greschner


% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('bin', 10);% ms
    p.addParamValue('corrcoef', true);% corrcoef or correlation index
    p.addParamValue('synchrony_index_field', 'synchrony_index');%
    p.addParamValue('verbose', true);% 
    % parse inputs
    p.parse(varargin{:});
    params = p.Results;
 
    
if params.verbose
    disp('activity movie - running');
end
  

%get synchrony_index        
if isfield(datarun, params.synchrony_index_field);
    synchrony_index=datarun.(params.synchrony_index_field);
else
    synchrony_index=zeros(length(datarun.cell_ids));
    synchrony_index=sparse(synchrony_index);
end
        
        
% get cell numbers
index_1=get_cell_indices(datarun,cell_specification_1);
if ~exist('cell_specification_2','var');
    index_2=index_1;
else
    index_2 = get_cell_indices(datarun,cell_specification_2);
end




%%%%corrcoef%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if params.corrcoef
     
    %split for speed vs memory 
    if isequal(index_1, index_2) 
        try
            spiketrain=zeros(datarun.duration*1000/params.bin,length(index_1));
            for i=1:length(index_1)
                spiketrain(:,i)=histc(datarun.spikes{index_1(i)},[0:params.bin/1000:datarun.duration-params.bin/1000]);    
            end
            temp=corrcoef(spiketrain);
            synchrony_index(index_1, index_2)=temp;
            fast=true;
        catch
            fast=false;
        end
    else
       try 
            spiketrain=zeros(datarun.duration*1000/params.bin,length(index_1)+length(index_2));
            for i=1:length(index_1)
                spiketrain(:,i)=histc(datarun.spikes{index_1(i)},[0:params.bin/1000:datarun.duration-params.bin/1000]);    
            end
            for i=1:length(index_2)
                spiketrain(:,i+length(index_1))=histc(datarun.spikes{index_2(i)},[0:params.bin/1000:datarun.duration-params.bin/1000]);    
            end
            temp=corrcoef(spiketrain);
            temp=temp(1:length(index_1),length(index_1)+1:end);
            synchrony_index(index_1, index_2)=temp;
            synchrony_index(index_2, index_1)=temp';
            fast=true;
       catch
            fast=false;
       end
    end
    
    if ~fast
        if params.verbose
            disp('memory limited - use loop');
            T=text_waitbar('correlation');
        end
        spiketrain=zeros(datarun.duration*1000/params.bin,2);
        for i=1:length(index_1)
            spiketrain(:,1)=histc(datarun.spikes{index_1(i)},[0:params.bin/1000:datarun.duration-params.bin/1000]); 
            for ii=1:length(index_2) 
                spiketrain(:,2)=histc(datarun.spikes{index_2(ii)},[0:params.bin/1000:datarun.duration-params.bin/1000]); 
                tem=corrcoef(spiketrain);
                temp(i,ii)=tem(1,2);
            end
            if params.verbose
                T=text_waitbar(T,i/length(index_1));
            end
        end
        synchrony_index(index_1, index_2)=temp;
        synchrony_index(index_2, index_1)=temp';
    end
       
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    

datarun.synchrony_index=synchrony_index;    
 
    
  
  
    
    
    


















    
    