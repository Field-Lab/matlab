function sp = get_spike_segment(spikes,tmin,tmax,type)


sp = spikes(spikes>tmin & spikes<tmax);

if (nargin > 3 && strcmp(type,'rel'))
    sp = sp - tmin;
end