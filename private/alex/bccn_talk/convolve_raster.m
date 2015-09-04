function convolve_raster(spikes, sig)



sr=120;
sig=50;
st=10000/sr*6.2/(60*sig);
time_list=-3.1:st:3.1;
kern=zeros(1,length(time_list));
for i=1:length(time_list)    
    kern(i)=250/sig*exp((1-time_list(i)^2)/2);
end
