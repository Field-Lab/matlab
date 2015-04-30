function [nspace ntime] = get_nspacetimeAH(Basepars)

switch(Basepars.k_filtermode)
                            case {'sep_raw','rk2'}
                                nspace = Basepars.k_spacepixels;
                                ntime = Basepars.k_stimframes;
                            case 'fixfilt' % added, edoi, 2012-01-04
                                nspace = Basepars.k_spacepixels;
                                ntime = Basepars.k_stimframes;
                            case 'sep_basis'
                                nspace = Basepars.nofilters_kspace;
                                ntime = Basepars.nofilters_ktime;
                            case {'nonsep','raw'}
                                error('Cannot get nspace or ntime in nonsep mode!');
                            end
