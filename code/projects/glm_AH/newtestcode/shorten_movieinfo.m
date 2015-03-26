
% hack generate shorter movie information files 


stimdir = '/Users/akheitman/NSEM_Home/Stimuli'
GLMType.fit_type = 'NSEM';


for i_type = 1:2
    if i_type == 1, GLMType.cone_model = '8pix_Identity_8pix';   end
    if i_type == 2, GLMType.cone_model = '8pix_Model1_1e4_8pix'; end
    
    for i_exp = 2:4
        if i_exp == 1, exp_nm = '2012-08-09-3'; end
        if i_exp == 2, exp_nm = '2012-09-27-3'; end
        if i_exp == 3, exp_nm = '2013-08-19-6'; end
        if i_exp == 4, exp_nm = '2013-10-10-0'; end
    
        [~, fullinputstats, ~] = loadmoviematfile(exp_nm , GLMType.fit_type, GLMType.cone_model,'fitmovie');
        
        a.stimulus       = fullinputstats.stimulus;
        a.cmodel         = fullinputstats.cmodel;
        a.dim            = fullinputstats.dim;
        a.mu_avgIperpix  = 255 * fullinputstats.mu_avgIperpix;
        %a.std_avgIperpix = 255 * fullinputstats.std_avgIperpix;
        a.range          = 255;
        a.hist_8bit      = fullinputstats.hist_8bit;
        a.histnote       = 'note the first element of the histogram is the value 0, 26th element is the value 255';
        
        clear inputstats
        inputstats = a;
        
        savedir = sprintf('%s/NSEM_%s', stimdir, a.stimulus);
        eval(sprintf('save %s/short_inputstats_%s.mat inputstats', savedir, inputstats.cmodel))
        
    end
end



%%%%
