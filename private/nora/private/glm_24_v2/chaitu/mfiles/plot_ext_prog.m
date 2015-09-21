function out = plot_ext_prog(x,j,basepars,stimpars,trainpars)

out = 0;

cla;


[taxis ext_x] = plot_ext_signal(x,basepars,stimpars.dt,size(trainpars.D,2));

sp_idx = (basepars.frame_offset)*(basepars.fac)+1:(basepars.frame_offset+basepars.maxt)*basepars.fac;
1;
[blah r] = plot_empirical_rate(trainpars.D,trainpars.dt,1,basepars.maxt*stimpars.dt,1);

r = r-mean(r);
r = r./std(r)*std(ext_x);

hold on, plot(makeaxis(1,length(r),(basepars.frame_offset-1)*stimpars.dt),r,'b');
xlim([basepars.frame_offset*stimpars.dt (basepars.frame_offset+basepars.maxt)*stimpars.dt]);