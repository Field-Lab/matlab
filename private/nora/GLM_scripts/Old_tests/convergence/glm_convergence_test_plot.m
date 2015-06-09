base_path='/Volumes/Analysis/nora/NSEM/GLM_Output/debug_fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams';
%fit_types={'WN','NSEM'};
%exp={'2012-08-09-3','2012-09-27-3','2013-08-19-6','2013-10-10-0'};
cell_count=0;
fit=dir(base_path);
for i=1:length(fit)
    if ~strcmp(fit(i).name(1),'.')
        expnm=dir([base_path '/' fit(i).name]);
        for j=1:length(expnm)
            if ~strcmp(expnm(j).name(1),'.')
                lengths=dir([base_path '/' fit(i).name '/' expnm(j).name]);
                for k=1:length(lengths)
                    if ~strcmp(lengths(k).name(1),'.')
                        cells=dir([base_path '/' fit(i).name '/' expnm(j).name '/' lengths(k).name]);
                        for l=1:length(cells)
                            if strcmp(cells(l).name(1), 'O')
                                cell_count=cell_count+1;
                                load([base_path '/' fit(i).name '/' expnm(j).name '/' lengths(k).name '/' cells(l).name]);
                                cellname=[ 'cell{' num2str(cell_count) '}.' fit(i).name];
                                eval([cellname '.BPS(' num2str(k) ')=fittedGLM.xvalperformance.glm_normedbits;'])
                                
                                debug_number=str2double(lengths(k).name((end-1):end));
                                    if isempty(debug_number)
                                        debug_number=str2double(lengths(k).name(end));
                                    end
                                    
                                eval([cellname '.debug=' num2str(debug_number)]);
                                disp([i j k l])
                            end
                        end
                    end
                end
            end
        end
    end
end

