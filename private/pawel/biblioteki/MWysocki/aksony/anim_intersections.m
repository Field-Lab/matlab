PatternNumber = 498; 
MovieNumber = 49:2:63;

figure(1)
filename = 'anim.gif';
for n = 1:length(MovieNumber)
      [wykresy,odl] = propagDirection(PatternNumber,MovieNumber(n));
      for i = 1:length(odl)
        subplot(2,3,i)
        plot(wykresy(i,:,1),wykresy(i,:,2)/max(wykresy(i,:,2)));
        xlabel('r');
        ylabel('amplitude');
        title(['Amplituda = ' num2str(MovieNumber(n)) ', odl = '  num2str(odl(i))]);
      end
      drawnow
      frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if n == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end
end
