%lines61 = load('output/m61_lines.mat');
lines63 = load('output/m63_lines.mat');
%lines59 = load('output/m59_lines.mat');
%lines = cat(1,lines61.lines,lines63.lines,lines59.lines);
lines = lines63.lines;

PatternNumber = [1:64,321:384]; 

x = -945:10:945;

electrodeMap = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);

figure(1)
%subplot(2,1,1)
for i = 1:length(lines(:,1))
        y = polyval(lines(i,:),x);
        plot(x,y);
        hold on
        axis([-945 945 -450 450])
end

for i = 1:length(lines(:,1))
        ex = electrodeMap.getXPosition(PatternNumber(i));
        ey = electrodeMap.getYPosition(PatternNumber(i));
        v = [100;polyval(lines(i,:),100)-polyval(lines(i,:),0)];
        v_norm = sqrt(v'*v);
        v = v/v_norm;
        przes = [-lines(i,2)/lines(i,1);0];
        v0 = v*(v'*([ex;ey] - przes)) + przes;
        
       h = plot([ex v0(1)],[ey v0(2)],'y')
       set(h,'LineWidth',4)
        hold on
        scatter(ex,ey,'r')
        axis([-945 945 -450 450])
end
% subplot(2,1,2)
% lines = lines63.lines;
% for i = 1:length(lines(:,1))
%         y = polyval(lines(i,:),x);
%         plot(x,y);
%         hold on
%         axis([-945 945 -450 450])
% end