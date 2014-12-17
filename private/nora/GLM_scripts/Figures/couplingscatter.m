clear
datapath='/Users/Nora/Desktop/research/GLM_Output';
CPpaths{1}='/fixedSP_rk1_linear_MU_PS_CP_p8IDp8/standardparams/';
CPpaths{2}='/fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
fittypepath{1}='NSEM_mapPRJ/';
fittypepath{2}='WN_mapPRJ/';

count=1;
for coupling=1:2
    for fittype=1:2
        for exp=1:4
            % filelist=ls([datapath CPpaths{coupling} fittypepath{fittype} exp_names(exp,:) '*.mat']);
            matfiles=dir([datapath CPpaths{coupling} fittypepath{fittype} exp_names(exp,:) '*.mat']);
            for file=1:length(matfiles)
                load([datapath CPpaths{coupling} fittypepath{fittype} exp_names(exp,:) matfiles(file).name]);
                bps(count)=fittedGLM.xvalperformance.glm_normedbits;
                count=count+1;
                % cell(coupling,fittype,exp,file)=filelist(matfiles(file).name);
            end
        end
    end
end

bps(bps<0)=0;
scatt=scatter(bps(79:156),bps(1:78));
hold on
plot ([0 1], [0 1],'black');
ylabel('Coupling')
xlabel('No Coupling')
axis([0 0.7 0 0.7]);
set(scatt,'Marker','.','MarkerEdgeColor','black');
