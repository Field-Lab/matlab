function [taxis ext_x] = plot_ext_signal(p,basepars,dt,Neff)


1;
ext = p(get_pars_idx(basepars,1,Neff,'ext'));


[taxis ext_x] = get_ext_signal(ext,basepars,dt);
1;
cla, plot(dt.*(basepars.frame_offset+ basepars.ext_timepts),ext,'ko','LineWidth',4);
hold on, plot(taxis,ext_x,'r','LineWidth',2);
title('Extrinsic signal'), xlabel('time (s)');
