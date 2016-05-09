function t_frame = t_frame_interpAH(t_trig,interp_method,fl_fig)
% return the timing of stimulus frames, based on 100-frames trigger.
% Usage: t_frame = t_frame_interp(t_trig,interp_method,fl_fig)

% AKHeitman  inheriting and modifying 2013-12-04
% edoi@salk.edu, 2011-11-16.

if ~exist('interp_method','var')
    interp_method = 'linear';
end
if ~exist('fl_fig','var')
    fl_fig = 0;
end

N = 100; % number of frames per trigger
n_frame = (length(t_trig)-1)*N+1; % total no of frames
idx_frame = 0:(n_frame-1);   % subtract 1 because the index starts with 0.
idx_trig  = 0:N:(n_frame-1); % ditto.
if idx_frame(end) ~= idx_trig(end), fprintf('error\n'), end

t_frame = interp1(idx_trig,t_trig,idx_frame,interp_method); 
t_frame = t_frame(1:(end-1));

% check with visualization 
if fl_fig
    figure(1), clf
    plot(idx_trig,t_trig,'r.');
    hold on
    plot(idx_frame(1:(end-1)),t_frame,'b')
    xlabel(sprintf('frame idx [0,%d]',idx_trig(end)))
    xlim([0,idx_trig(end)])
    ylabel('time [sec]')
    title('Frame time')
    
    figure(2), clf
    plot(idx_frame(1:(end-2)),diff(t_frame),'b.')
    xlim([0,idx_trig(end-1)])
    xlabel(sprintf('frame idx [0,%d]',idx_trig(end)-1))
    ylabel('time [sec]')
    title('Difference of frame time')

    figure(3), clf
    hist(diff(t_trig))
    title('Histogram of difference of triggers')
    % the jitter is likely to be in the unit of sampling interval, i.e.,
    % 1/datarun.sampling_rate.
end
