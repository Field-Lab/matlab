

function [X stim_params] = get_wn_stimulus(datafile,totalframes)

if (strfind(datafile,'2005'))
    stixel_val = 20;
    stim_params = struct('mdf','/e/1.3/p1/chaitu/data/rrs-library-3.1/mdf/BW-20-8-0.48-11111.mdf', ...
        'color','grey','stixel',stixel_val);
elseif (strfind(datafile,'2008-06'))
    %BW-10-16-0.48,
    stixel_val = 10;
    stim_params = struct('mdf','/e/1.3/p1/chaitu/data/rrs-library-3.1/mdf/BW-10-16-0.48-11111.mdf', ...
        'color','grey','stixel',stixel_val);
elseif(strfind(datafile,'2008-08'))
    %BW-10-8-0.48
    stixel_val = 10;
    stim_params = struct('mdf','/e/1.3/p1/chaitu/data/rrs-library-3.1/mdf/BW-10-8-0.48-11111.mdf', ...
        'color','grey','stixel',stixel_val);
elseif(strfind(datafile,'2009-11'))
    1;
    stixel_val = 8;
    stim_params = struct('mdf','/e/1.3/p1/chaitu/java/mdf/BW-8-8-0.48-11111.mdf',...
        'color','grey','stixel',stixel_val);
elseif(strfind(datafile,'2009-12'))
    1;
    stixel_val = 8;
    stim_params = struct('mdf','/e/1.3/p1/chaitu/java/mdf/BW-16-8-0.48-11111.mdf',...
        'color','grey','stixel',stixel_val);
elseif(strfind(datafile,'2010-03'))
    stixel_val = 8;
    stim_params = struct('mdf','/e/1.3/p1/chaitu/java/mdf/BW-8-8-0.48-11111.mdf',...
        'color','grey','stixel',stixel_val);
end
    
1;
addpath ~/data/rrs-library-3.1/
javaaddpath('/e/1.3/p1/chaitu/data/rrs-library-3.1/rrs.jar');
javaaddpath('/e/1.3/p1/chaitu/data/rrs-library-3.1/cell-finder.jar');
X = export_stimulus('wn', totalframes, [], stim_params)';