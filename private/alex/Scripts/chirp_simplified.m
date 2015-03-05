clear

chirpLength=22500; % length of the chirp stimulus you want to take
date='20121023';

% formatted chirp data
chirpPath=dir(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\*chirp*']);
load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',chirpPath(1).name])
trialsChirp=size(spike_info.name_info,1);

for unit=1:length(chirpPath)
    load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',name,'___spike_info_chirp'])
    
    chirp=zeros(trialsChirp,25000);
    for trial=1:trialsChirp
        spikes=spike_info.spike_times{trial};
        flips=spike_info.flip_times{trial}(end,1);
        conv=convolved(spikes,40,flips);
        conv=conv(121:end-120);
        chirp(trial,1:length(conv))=conv;
    end
    cnt=1;
    % I had 5 repetitions of chirp per ND, so I average every 5 responses
    for i=1:5:35
        responses(cnt,1:chirpLength,unit)=mean(chirp(i:i+4,1:chirpLength));
        cnt=cnt+1;
    end
end

%% SLIDING WINDOW

noStim=200:2500;
window=500;
% structure: rows - parameters,columns - nds, slices - units
for i=1:35
    for j=1:7
        baseLine=mean(responses(j,noStim,i));        
        based=abs(responses(j,:,i)-baseLine);
% based=responses(j,:,i);
        spont(j,i)=sum(based(noStim));
        cnt=1;
        for fc=1:50:chirpLength-window
            param(cnt,j,i)=sum(based(fc:fc+window))/window;
            cnt=cnt+1;
        end        
        spont(1:cnt-1,j,i)=repmat(sum(based(noStim)),cnt-1,1);
    end
end

%param is the output variable - "modulation strength"