function b1 = loop_slider(k,kmin,kmax)

   %b1 = loop_slider(1,kmin,kmax);
   %k=get(b1,'Value') 
     
   b1= uicontrol(gcf,...
      'Style','slider',...
      'Min' ,kmin,'Max',kmax, ...
      'Position',[0,0,200,15], ...
      'Value', k,...
      'SliderStep',[1/(kmax-kmin) 1/(kmax-kmin)],...
      'CallBack', 'uiresume;');
    
  
 