% RASTERPLOT.M Display spike rasters.
%   RASTERPLOT(T,N,L) Plots the rasters of spiketimes (T in samples) for N trials, each of length
%   L samples, Sampling rate = 1kHz. Spiketimes are hashed by the trial length.
% 
%   RASTERPLOT(T,N,L,H) Plots the rasters in the axis handle H
%
%   RASTERPLOT(T,N,L,H,FS) Plots the rasters in the axis handle H. Uses sampling rate of FS (Hz)
%
%   Example:
%          t=[10 250 9000 1300,1600,2405,2900];
%          rasterplot(t,3,1000)
%

% Rajiv Narayan
% askrajiv@gmail.com
% Boston University, Boston, MA

function rasterplot_options(times,numtrials,triallen, starter, plotcolor, hresp)


%%%%%%%%%%%%%% Plot variables %%%%%%%%%%%%%%
plotwidth=1;     % spike thickness
trialgap=1.5;    % distance between trials
defaultfs=1000;  % default sampling rate
showtimescale=1; % display timescale
showlabels=0;    % display x and y labels

fs=defaultfs;

 % plot spikes
 
 trials=ceil(times/triallen);
 reltimes=mod(times,triallen);
 reltimes(~reltimes)=triallen;
 numspikes=length(times);
 xx=ones(3*numspikes,1)*nan;
 yy=ones(3*numspikes,1)*nan;
 
 yy(1:3:3*numspikes)=(trials-1)*trialgap;
 yy(2:3:3*numspikes)=yy(1:3:3*numspikes)+1;
 yy = yy + starter;
 
 %scale the time axis to ms
 xx(1:3:3*numspikes)=reltimes*1000/fs;
 xx(2:3:3*numspikes)=reltimes*1000/fs;
 xlim=[1,triallen*1000/fs];
 
 axes(hresp);
 plot(xx, yy, plotcolor, 'linewidth',plotwidth);
 axis ([xlim,0,(numtrials)*1.5]);
 
 if (showtimescale)
     set(hresp, 'ytick', [],'tickdir','out');
 else
     set(hresp,'ytick',[],'xtick',[]);
 end
 
 if (showlabels)
     xlabel('Time(ms)');
     ylabel('Trials');
 end