function spks=binSpks(spikes,fr)


%spikes=datarun.spikes{cellID};
% spikes are in s. convert to ms
spikes=round(spikes*1000);
spks=zeros(length(fr)-1,1);
for ibin=1:length(fr)-1
spks(ibin)=sum(double(spikes>=fr(ibin) & spikes<fr(ibin+1)));
end

end