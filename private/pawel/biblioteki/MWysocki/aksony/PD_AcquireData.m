function [ file ] = PD_AcquireData(InputPath, OutputPath, DataID, PatternRange, MovieRange, BadElectrodes, InitialDirection)

    file = cell(length(MovieRange));
    prefile = [OutputPath filesep 'PD_' DataID '_' datestr(now,'yyyy_mm_dd_HH_MM') '_m'];
    
    line = cell(length(PatternRange),1);
    orient = cell(length(PatternRange),1);
    v0 = cell(length(PatternRange),1);
    ve = cell(length(PatternRange),1);
    vp = cell(length(PatternRange),1);
    edges = cell(length(PatternRange),1);
    func = cell(length(PatternRange),1);
    success = cell(length(PatternRange),1);
    rms = cell(length(PatternRange),1);
    
    patterns = PatternRange;
    
    for m = 1:length(MovieRange)
        for p = 1:length(PatternRange)
            [line0, orient0, v00, ve0, vp0, edges0, func0, success0, rms0] = PD_GetPropagDirection (InputPath,PatternRange(p),MovieRange(m),BadElectrodes,InitialDirection);
            line{p} = line0;
            orient{p} = orient0;
            v0{p} = v00;
            ve{p} = ve0;
            vp{p} = vp0;
            edges{p} = edges0;
            func{p} = func0;
            success{p} = success0;
            rms{p} = rms0;
        end
        file{m} = [prefile num2str(MovieRange(m)) '.mat'];
        save(file{m},'line','orient','v0','ve','vp','edges','func','success','rms','patterns');
    end
end

