function cellInfo = cell_list_sbc(exclude30ArrayData)

ii = 0;

%complete analysis (responses analyzable through p>0.5)
ii = ii + 1;
cellInfo(ii).id             = 466; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 28;
cellInfo(ii).verMin         = true; %true IFF stim elec not on edge of array or next to dead elec (4) AND all surrounding stim elecs have verified higher thresholds
cellInfo(ii).analysisPath   = '2008-08-26-0/data002';
cellInfo(ii).PW             = 50;
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 148.9; %in square microns

ii = ii + 1;
cellInfo(ii).id             = 662; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 49;
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2010-03-05-0/data001';
cellInfo(ii).PW             = 100;
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 194.9;

ii = ii + 1;
cellInfo(ii).id             = 873; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 59;
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2010-10-18-3/data010';
cellInfo(ii).PW             = 100;
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 210.9;

ii = ii + 1;
cellInfo(ii).id             = 18; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 5;
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2010-10-28-2/data003';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 197.6;

ii = ii + 1;
cellInfo(ii).id             = 242; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 15;
cellInfo(ii).verMin         = true;
cellInfo(ii).analysisPath   = '2010-10-28-3/data002';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 30;
cellInfo(ii).meanElecArea   = [];

ii = ii + 1;
cellInfo(ii).id             = 544; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 37;
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2010-10-28-3/data002';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 30;
cellInfo(ii).meanElecArea   = [];

ii = ii + 1;
cellInfo(ii).id             = 376; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 23;
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2010-11-22-3/data002';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 30;
cellInfo(ii).meanElecArea   = [];

ii = ii + 1;
cellInfo(ii).id             = 586; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 40;
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2010-11-22-3/data002';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 30;
cellInfo(ii).meanElecArea   = [];

ii = ii + 1;
cellInfo(ii).id             = 902; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 61;
cellInfo(ii).verMin         = true;
cellInfo(ii).analysisPath   = '2011-01-11-0/data002_003';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 166.3;

ii = ii + 1;
cellInfo(ii).id             = 798; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 47;
cellInfo(ii).verMin         = true;
cellInfo(ii).analysisPath   = '2011-05-11-2/data002';
cellInfo(ii).PW             = 100;
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 200.7;

ii = ii + 1;
cellInfo(ii).id             = 346; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 24;
cellInfo(ii).verMin         = true;
cellInfo(ii).analysisPath   = '2011-05-11-5/data002';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 30;
cellInfo(ii).meanElecArea   = [];

ii = ii + 1;
cellInfo(ii).id             = 768; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 52;
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2011-06-24-0/data003';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 116.1;

ii = ii + 1;
cellInfo(ii).id             = 271; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 19;
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2011-07-05-0/data002';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 115.5;

ii = ii + 1;
cellInfo(ii).id             = 346; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 29;
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2011-07-05-0/data002';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 115.5;

ii = ii + 1;
cellInfo(ii).id             = 662; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 45;
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2011-07-05-0/data002';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 115.5;

ii = ii + 1;
cellInfo(ii).id             = 393; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 27;
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2011-08-04-0/data002_003';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 126.4;

ii = ii + 1;
cellInfo(ii).id             = 136; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 7;
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2011-10-25-4/data002';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 128.7;

ii = ii + 1;
cellInfo(ii).id             = 334; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 21;
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2011-10-25-4/data002';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 128.7;

ii = ii + 1;
cellInfo(ii).id             = 632; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 46;
cellInfo(ii).verMin         = true;
cellInfo(ii).analysisPath   = '2011-10-25-4/data002';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 128.7;

ii = ii + 1;
cellInfo(ii).id             = 93; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 5;
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2012-01-27-3/data002';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 135.0;

ii = ii + 1;
cellInfo(ii).id             = 241; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 12;
cellInfo(ii).verMin         = true;
cellInfo(ii).analysisPath   = '2012-01-27-3/data002';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 135.0;

ii = ii + 1;
cellInfo(ii).id             = 331; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 20;
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2012-01-27-3/data002';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 135.0;

ii = ii + 1;
cellInfo(ii).id             = 423; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 29;
cellInfo(ii).verMin         = true;
cellInfo(ii).analysisPath   = '2012-01-27-3/data002';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 135.0;

%%% partial analysis (responses but not analyzable through p>0.5)

ii = ii + 1;
cellInfo(ii).id             = 796; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 55;
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2008-11-12-0/data001';
cellInfo(ii).PW             = 50;
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 30;
cellInfo(ii).meanElecArea   = [];

ii = ii + 1;
cellInfo(ii).id             = 197; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 12; %only analyzable through p = 0.42 (analyzable for more of curve for stim elec 6 but has higher thresh)
cellInfo(ii).verMin         = true;
cellInfo(ii).analysisPath   = '2010-10-18-3/data001';
cellInfo(ii).PW             = 100;
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 210.9;

ii = ii + 1;
cellInfo(ii).id             = 481; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 33;
cellInfo(ii).verMin         = true;
cellInfo(ii).analysisPath   = '2010-10-28-2/data003';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = 197.6;

ii = ii + 1;
cellInfo(ii).id             = 932; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 63;
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2011-05-11-5/data002';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 30;
cellInfo(ii).meanElecArea   = [];

%incomplete analysis / no stimulation detected
% min thresh values are based on fitting an erf to analyzed movies +
% hypothetical p = 1 movies at 10% increase amplitude increments

ii = ii + 1;
cellInfo(ii).id             = 572; %#ok<*SAGROW>
cellInfo(ii).stimElec       = [];
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2008-11-10-3/data001';
cellInfo(ii).limStimElec    = 39;
cellInfo(ii).PW             = 50;
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = [];

ii = ii + 1;
cellInfo(ii).id             = 661; %#ok<*SAGROW>
cellInfo(ii).stimElec       = [];
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2008-11-10-3/data001';
cellInfo(ii).limStimElec    = 45;
cellInfo(ii).PW             = 50;
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = [];

ii = ii + 1;
cellInfo(ii).id             = 873; %#ok<*SAGROW>
cellInfo(ii).stimElec       = [];
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2008-11-12-3/data001';
cellInfo(ii).limStimElec    = 61;
cellInfo(ii).PW             = 50;
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 30;
cellInfo(ii).meanElecArea   = [];

ii = ii + 1;
cellInfo(ii).id             = 557; %#ok<*SAGROW>
cellInfo(ii).stimElec       = [];
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2011-10-25-3/data001';
cellInfo(ii).limStimElec    = 33;
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 30;
cellInfo(ii).meanElecArea   = [];

ii = ii + 1;
cellInfo(ii).id             = 769; %#ok<*SAGROW>
cellInfo(ii).stimElec       = [];
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2011-10-25-4/data002';
cellInfo(ii).limStimElec    = 51;
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 60;
cellInfo(ii).meanElecArea   = [];

% haven't attempted to analyze yet because probably doesn't meet ei
% amplitude threshold or is on 30 µm array
ii = ii + 1;
cellInfo(ii).id             = 766; %#ok<*SAGROW>
cellInfo(ii).stimElec       = 0; %indicates unfinished analysis
cellInfo(ii).verMin         = false;
cellInfo(ii).analysisPath   = '2011-10-25-3/data001';
cellInfo(ii).PW             = [50 100];
cellInfo(ii).type           = 'sbc';
cellInfo(ii).arraySpacing   = 30;
cellInfo(ii).meanElecArea   = [];

if exclude30ArrayData
    % remove cells that were on 30 micron arrays
    for ii = length(cellInfo):-1:1
        if cellInfo(ii).arraySpacing == 30;
            cellInfo(ii) = [];
        elseif cellInfo(ii).arraySpacing ~= 60;
            error(['arraySpacing should be 30 or 60 but it''s' num2str(cellInfo(ii).arraySpacing)])
        end
    end
end



