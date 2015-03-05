function convolve=convolved(spikes,sig,tt)
% returns firing rate convolving a spike sequence with a gaussian. 
% spikes - sequence of spike times (ms), sampling frequency then 1kHz
% sig - sigma of the gaussian kernel
% tt - fill with 0 till time = tt. If tt=0, last spike time is the last
% time point
if tt==0
    tt=spikes(end);
end
sr=1000;
st=10000/sr*6.2/(60*sig);
time_list=-3.1:st:3.1;
kern=zeros(1,length(time_list));
for i=1:length(time_list)    
    kern(i)=250/sig*exp((1-time_list(i)^2)/2);
end

channel=ceil(spikes);
channel(channel>tt)=[];
binary_spikes=zeros(1,tt);
channel(channel<=0)=[];
binary_spikes(channel)=1;
binary_spikes(channel(diff(channel)==0))=2;
convolve=conv(binary_spikes,kern);

end