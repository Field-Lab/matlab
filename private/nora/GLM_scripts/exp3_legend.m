function exp3_legend(axes_handle, short)
if short
    strings = {'On', 'Off'};
else
    strings = {'2012-08-09-3 ON','2012-08-09-3 OFF', ...
        '2012-09-27-3 ON','2012-09-27-3 OFF', ...
        '2013-08-19-6 ON', '2013-08-19-6 OFF'};
end
legend(axes_handle, strings)
end