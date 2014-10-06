%h=NS512_ShowEIFrameAsCircles2(EI,ArrayID,Channels,MarkedChannels,Colors1,MarkedChannels2,Colors2,XRange,YRange);
EI=rand(512,100);
for i=1:100
    h=NS512_ShowEIAsCircles(EI(,i),500,[1:512],[],[-1000 1000],[-500 500]);
    pause(0.2)
end
