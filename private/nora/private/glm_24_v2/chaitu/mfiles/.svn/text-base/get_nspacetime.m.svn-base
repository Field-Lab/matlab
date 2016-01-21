function [nspace ntime] = get_nspacetime(basepars)

switch(basepars.filtermode)
    case {'sep_raw','rk2'}
        nspace = basepars.n;
        ntime = basepars.Mk;
    case 'fixfilt' % added, edoi, 2012-01-04
        nspace = basepars.n;
        ntime = basepars.Mk;
    case 'sep_basis'
        nspace = basepars.nofilters_kspace;
        ntime = basepars.nofilters_ktime;
    case {'nonsep','raw'}
        error('Cannot get nspace or ntime in nonsep mode!');
end
