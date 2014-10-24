function [I_sc dat] = mkSubAssgn(C,type)

%-------------------------------------------
%
% out = mkSubAssgn(sim)
%
% create binary subunit assignments for 
% simulating or fitting an LN cascade model
% 
% inputs
% C -- number of stimuli
% type -- form of assignments
%         options are:
%           'singletons'
%           'doublets'
%           'doublets-overlap'
%           'singletons-and-doublets'
%
% outputs
% I_sc -- desired subunit assignment matrix
% dat -- dat structure with fileds filled in
% 
% freeman, 3-19-2012
%-------------------------------------------

% location of saved structures for simulation
path = '~/Dropbox/Projects/retinaSubunits/examples';

% set cone locations to empty
% (only used for some structures)
dat.locs_c = [];

switch type
    case 'singletons'
        % all singletons
        I_sc = eye(C);
    case 'doublets'
        % all doublets, return an error if not divisible by 2
        if mod(C,2) ~= 0
            error('(mkSimSubunits) number of stimuli must divide by 2 to make doublets');
        end
        S = C/2;
        I_sc = zeros(S,C);
        for is=1:S
            I_sc(is,(is*2-1):is*2) = 1;
        end
    case 'doublets-overlap'
        % all doublets, return an error if not divisible by 2
        S = C - 1;
        I_sc = zeros(S,C);
        for is=1:S
            I_sc(is,is:is+1) = 1;
        end
    case 'singletons-and-doublets'
        % half doublets, half singletons
        nDouble = S/2;
        nSingle = S/2;
        if (nDouble*2 + nSingle) ~= C
            error('(mkSimSubunits) mismatch between number of stimuli and number of subunits');
        end
        I_sc = zeros(S,C);
        for is=1:nDouble
            I_sc(is,(is*2-1):is*2) = 1;
        end
        for is = nDouble+1:nDouble+nSingle
            I_sc(is,is+nDouble) = 1;
        end
    case 'midget'
        load(fullfile(path,'midget-sim.mat'));  
        I_sc = fit.I_sc;
        dat.locs_c = fit.locs_c;
    case 'midget-overlap'
        load(fullfile(path,'midget-overlap-sim.mat'));  
        I_sc = fit.I_sc;
        dat.locs_c = fit.locs_c;
    case 'midget-singletons'
        load(fullfile(path,'midget-overlap-sim.mat'));  
        I_sc = eye(size(fit.locs_c,1),size(fit.locs_c,1));
        dat.locs_c = fit.locs_c;
    case 'parasol'
        load(fullfile(path,'parasol-sim.mat'));  
        I_sc = fit.I_sc;
        dat.locs_c = fit.locs_c;
    case 'parasol-overlap'
        load(fullfile(path,'parasol-overlap-sim.mat'));  
        I_sc = fit.I_sc;
        dat.locs_c = fit.locs_c;
    case 'parasol-singletons'
        load(fullfile(path,'parasol-overlap-sim.mat')); 
        I_sc = eye(size(fit.locs_c,1),size(fit.locs_c,1));
        dat.locs_c = fit.locs_c;
end

if ~isempty(dat.locs_c)
    dat.xRange = round([min(dat.locs_c(:,1))*0.98 max(dat.locs_c(:,1))*1.02]);
    dat.yRange = round([min(dat.locs_c(:,2))*0.98 max(dat.locs_c(:,2))*1.02]);
end


dat.cellType = 'simulation';
dat.dataSet = 'simulation';
dat.rgcId = 0;