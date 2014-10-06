 function p=NS512_VanishingDataBySigmaKB(t,p,cf);
    switch logical(true)
            
        case (cf.sigma<1)
            
          t_min=max(round(cf.tau-(20*cf.sigma)),1);
          t_max=min(round(cf.tau+(20*cf.sigma)), t(length(t)));
        
           p(t_min:t_max)=0;
             
        case cf.sigma<2  
            
          t_min=max(round(cf.tau-(12*cf.sigma)),1);
          t_max=min(round(cf.tau+(12*cf.sigma)),  t(length(t)));
        
           p(t_min:t_max)=0;
         
        case cf.sigma<10
            
         t_min=max(round(cf.tau-(10*cf.sigma)),1);
         t_max=min(round(cf.tau+(10*cf.sigma)),  t(length(t)));
        
            p(t_min:t_max)=0;
        
        case cf.sigma <50
            
         t_min=max(round(cf.tau-(6*cf.sigma)),1);
         t_max=min(round(cf.tau+(6*cf.sigma)), t(length(t)));
        
            p(t_min:t_max)=0;
        
        case cf.sigma <100 
            
         t_min=max(round(cf.tau-(4*cf.sigma)),1);
         t_max=min(round(cf.tau+(4*cf.sigma)), t(length(t)));
        
          p(t_min:t_max)=0;
        
        case cf.sigma >=100 
            
         t_min=max(round(cf.tau-(3*cf.sigma)),1);
         t_max=min(round(cf.tau+(3*cf.sigma)), t(length(t)));
        
         p(t_min:t_max)=0;
      
        otherwise warning('Unexpected sigma');

    end