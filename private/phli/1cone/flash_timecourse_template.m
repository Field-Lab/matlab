function [flashtemplate, flashtemplatex] = flash_timecourse_template(template, templaterun, flashstim)
% FLASH_TIMECOURSE_TEMPLATE    Convolve STA timecourse with a flash timecourse, handle time points
%
% 2012-06 phli, abstracted from 2012-04-13-1_analysis script
%

% Build timecourse response templates from earlier white noise run
% x scaled to fundamental monitor rate
interval = templaterun.stimulus.interval;
monitor_refresh = templaterun.stimulus.refresh_period / interval / 1000;
sta_offset = templaterun.stas.java_sta.getSTAOffset;
timecourse_template  = interp1(1:length(template), template, 1:1/interval:length(template), 'cubic', 0);
timecourse_templatex = ((1:1/interval:length(template)) - sta_offset - 1) .* monitor_refresh * interval;

% Convolve with length of flash for flash trial
flashcourse  = [0 ones(1, flashstim.spec.frames) 0];
flashcoursex = [-1:flashstim.spec.frames] .* monitor_refresh;
flashtemplate  = conv( timecourse_template,  flashcourse);
flashtemplatex = convx(timecourse_templatex, flashcoursex);