function [taxis ext_x] = get_ext_signal(ext,basepars,dt)


1;
taxis = makeaxis(dt,basepars.maxt,(basepars.frame_offset-1)*dt)';
if (~isfield(basepars,'interpfilt'))
    % Use sinc interpolation
    ext_x = sincinterp(makeaxis(dt,basepars.maxt),basepars.ext_timepts,ext);
else
    ext_x = geninterp(basepars.interpfilt,basepars.maxt,basepars.ext_timepts,ext);%sincinterp(taxis,basepars.ext_timepts, ext);
end
