% started AKHeitman2014-03-08  
% if we need cleaner data .. come here!


% 2014-04-26
% This is proviing to be difficult.
% Conclude .. doing by cel type across cells is impossible !! !!
% Maybe just do on a cell by cell basis.. some kind of algorithm?
% could give overall better results to knock out extra trials ..?


if strcmp(exp_nm, '2012-08-09-3') && (strcmp(stimtype, 'BW')  || strcmp(stimtype,'WN')  
    badrast_trials = [];
    OFFP_badrast_trials = [];
    ONP_badrast_trials = [];
end

if strcmp(exp_nm, '2012-08-09-3') && strcmp(stimtype, 'NSEM')
    badrast_trials = [];
    OFFP_badrast_trials = [];
    ONP_badrast_trials = [];
end


if strcmp(exp_nm, '2012-09-27-3') && (strcmp(stimtype, 'BW')  || strcmp(stimtype,'WN')  
    badrast_trials = [];
    OFFP_badrast_trials = [];
    ONP_badrast_trials = [25:32];
end