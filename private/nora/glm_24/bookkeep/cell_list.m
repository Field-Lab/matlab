% AKHeitman 2014-03-29
% Trying to bring better organization to our cell selection process
% Avoid annoying loop structures over and over again


function [exp_nm,cid_cell,expname]  = cell_list( expnumber, cellselectiontype)


%%% Identify exp_nm  %%%
if expnumber == 1, exp_nm ='2012-08-09-3'; expname = 'expA';  end
if expnumber == 2, exp_nm ='2012-09-27-3'; expname = 'expB';  end
if expnumber == 3, exp_nm ='2013-08-19-6'; expname = 'expC';  end 
if expnumber == 4, exp_nm ='2013-10-10-0'; expname = 'expD';  end 
if expnumber == 5, exp_nm ='2014-06-04-0'; expname = 'expE';  end 


%%% Find appropriate CID %%%
if strcmp(cellselectiontype, 'shortlist')    
    if expnumber == 1, cid_cell = {1471 3676  5086 5161 841 1426 1772 2101 1276}; end % 3676 and 1772 doesnt work for power raise on Bertha
    if expnumber == 2, cid_cell = {1 31 301 1201 1726 91 1909 2360 6858};  end
    if expnumber == 3, cid_cell = {2824 3167 3996 5660 6799 737 1328 1341 2959  5447}; end 
    if expnumber == 4, cid_cell =  {32 768 2778 4354 5866 7036 346 1233 3137  5042 5418};end %cid_cell = { 3137  5042 5418}
    if expnumber == 5, cid_cell = {901 5086}; end
end


if strcmp(cellselectiontype, 'debug')    
    if expnumber == 1, cid_cell = {841}; end %{5086 841 3676}
    if expnumber == 2, cid_cell = {1 6858};  end
    if expnumber == 3, cid_cell = {2824  5447}; end 
    if expnumber == 4, cid_cell = { 5866 5418};end
    if expnumber == 5, cid_cell = {901 5086}; end
end

if strcmp(cellselectiontype,'all')
    load('all_cells.mat');
    cid_cell = all_cells{expnumber};
end


%{
if ~exist('shead_cellID','var')
    if strcmp(GLMPars.fit_type , 'BW')
    head_cellID = [datarun_slv.cell_types{2}.cell_ids datarun_slv.cell_types{1}.cell_ids];% datarun_slv.cell_types{5}.cell_ids];
    elseif strcmp(GLMPars.fit_type , 'NSEM')
         [StimulusPars DirPars datarun_slvBW datarun_mas] = Directories_Params_v20_split(exp_nm, 'BW', map_type)
         head_cellIDBW = [datarun_slvBW.cell_types{1}.cell_ids datarun_slvBW.cell_types{2}.cell_ids];     
         head_cellID = [intersect(datarun_slv.cell_types{2}.cell_ids , head_cellIDBW) , intersect(datarun_slv.cell_types{1}.cell_ids , head_cellIDBW) ]
    end
end
%} 

end

%{
??? Error using ==> sfminbx at 206
fminunc cannot continue: user function is returning Inf or NaN values.

Error in ==> fminunc at 370
   [x,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN] = sfminbx(funfcn,x,l,u, ...

Error in ==> glm_DOUBLEOPT_innerloop at 66
        [pstar fstar eflag output]     = fminunc(@(p)
        glm_convex_optimizationfunction(p,glm_covariate_vec,home_spbins,t_bin),p_init,optim_struct);
        
Error in ==>
glm_DOUBLEOPT_execute>@(nonlin_par)glm_DOUBLEOPT_innerloop(nonlin_par,'input_pt_nonlinearity',glm_cellinfo,GLMType,GLMPars,fixed_covariates,home_spbins,fitmovie,t_bin,p_guess)
at 116
    [nonlinstar, fstar eflag output] = fminunc(@(nonlin_par)
    glm_DOUBLEOPT_innerloop(nonlin_par,'input_pt_nonlinearity',...

Error in ==> lineSearch>bracketingPhase at 103
          f_alpha =
          feval(funfcn{3},reshape(xInitial(:)+alpha*dir(:),sizes.xRows,sizes.xCols),varargin{:});
          
Error in ==> lineSearch at 48
[a,b,f_a,fPrime_a,f_b,fPrime_b,alpha,f_alpha,grad,exitflagBrckt,funcCountBrckt] = ...

Error in ==> fminusub at 208
    [alpha,f,grad,exitflagLnSrch,funcCountLnSrch] = ...

Error in ==> fminunc at 367
   [x,FVAL,GRAD,HESSIAN,EXITFLAG,OUTPUT] = fminusub(funfcn,x, ...

Error in ==> glm_DOUBLEOPT_execute at 116
    [nonlinstar, fstar eflag output] = fminunc(@(nonlin_par)
    glm_DOUBLEOPT_innerloop(nonlin_par,'input_pt_nonlinearity',...

Error in ==> glmwrap24_func at 142
                [fittedGLM] =    glm_DOUBLEOPT_execute(GLMType, spikesconcat,
                concat_fitmovie, glm_cellinfo);
 %}