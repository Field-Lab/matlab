function msg = genNeuronAnim
global gif_path gif_width gif_height gif_time_step param_path neuron_path neuron_classID1 neuron_classID2 double_mode time_res time_offset time_dur time_light;

%if length(input) ~= 8 | length(output) ~= 5
%    msg = {'Arguments lists are not valid! Aborting generation!'};
%    return;
%end

msg = {'Error...'};

%gif_path = output{1}
%gif_width = output{2};
%gif_height = output{3};
%gif_time_step = output{4};
%param_path = input{1};
%neuron_path = input{2};
%neuron_classID1 = input{3};
%neuron_classID2 = input{4};
%double_mode = input{5};
%time_res = input{6}; % [ms] (real time between frames)
%time_offset = input{7}; % [ms]
%time_dur = input{8}; % [ms]
%time_light = output{5}; %number of frames when neuron light on

time_offset_sec = time_offset/1000; % [s]

w = waitbar(0,'Generating...');

paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile(param_path);
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuron_path);
idList = paramsFile.getNeuronsInClass(neuron_classID1);
if double_mode
    idList2 = paramsFile.getNeuronsInClass(neuron_classID2);
end
neuron_freq = neuronFile.getSamplingFrequency();
%neuron_samples = neuronFile.getNumberOfSamples();
neuron_period = 1000/neuron_freq; % [ms]
%neuron_total_time = neuron_period*neuron_samples; % [ms]
samples_skip = time_offset/neuron_period;
samples_show = time_dur/neuron_period;
samples_frame = neuron_period/time_res; % 1/(number of samples in one frame)

t = 0:0.01:1; % parameter to draw ellipses
fig=figure;
set(fig,'Visible','off','color','w','Position',[0,0,gif_width,gif_height]);
%full_axes = axes('Position',[0 0 1 1],'Visible','off');
%graph_axes = axes('Box','on','Position',[.05 .1 .75 .85]); - for visible Ticks
if double_mode
    graph2_axes = subplot('Position',[0.42 0.04 0.39 0.91],'Box','on');
end
if double_mode
    graph_axes = subplot('Position',[0.02 0.04 0.39 0.91],'Box','on');
else
    graph_axes = subplot('Position',[0.02 0.04 0.78 0.91],'Box','on');
end
full_axes = axes('Position',[0 0 1 1],'Visible','off');
title(graph_axes,neuron_classID1);
set(graph_axes,'XTick',[]);
set(graph_axes,'YTick',[]);
if double_mode
    title(graph2_axes,neuron_classID2);
    set(graph2_axes,'XTick',[]);
    set(graph2_axes,'YTick',[]);
end

color = [[0 1 0];[1 0 0];[0 1 1]];
alpha = [0.95 0.95 0.95];

subplot(full_axes);
hold on
label = text(0.83,0.9,'Time:','FontSize',12);
time_text = text(0.9,0.9,[num2str(time_offset_sec,'%1.3f') ' s'],'FontSize',12,'BackgroundColor','w');
label2 = text(0.83,0.85,['Resolution:  ' num2str(time_res,'%1.2f') ' ms / frame'],'FontSize',12);
time_text = text(0.9,0.9,[num2str(time_offset_sec,'%1.3f') ' s'],'FontSize',12,'BackgroundColor','w');
hold off

subplot(graph_axes);
hold on
%M: nr_kola x0 y0 a b kat
%N: czas nr_neuronu kolor czas_zycia (czas zyci - na jak dlugo zapalic elipse)

M=zeros(length(idList),6);
N_tmp=[];
elli = zeros(length(idList),1);
for i=1:length(idList)
    neuronID=idList(i);
    M(i,1) = i;
    M(i,2) = paramsFile.getDoubleCell(neuronID(1), 'x0');
    M(i,3) = paramsFile.getDoubleCell(neuronID(1), 'y0');
    M(i,4) = paramsFile.getDoubleCell(neuronID(1), 'SigmaX');
    M(i,5) = paramsFile.getDoubleCell(neuronID(1), 'SigmaY');
    M(i,6) = paramsFile.getDoubleCell(neuronID(1), 'Theta');
    
    x1 = M(i,4)*sin(2*pi*t);
    y1 = M(i,5)*cos(2*pi*t);
    x2 = x1*cos(M(i,6)) + y1*sin(M(i,6)) + M(i,2);
    y2 = - x1*sin(M(i,6)) + y1*cos(M(i,6)) + M(i,3);
    if double_mode
       elli(i) = fill(y2,x2,'w','LineWidth',1);
    else
       elli(i) = fill(x2,y2,'w','LineWidth',1);
    end
    
    if neuronFile.containsID(neuronID)
        spikeTimes = neuronFile.getSpikeTimes(neuronID)';
        indeksy = find(spikeTimes > samples_skip & spikeTimes < (samples_skip + samples_show)); %skip time_offset ms and visualize time_dur ms
        czasy=round((spikeTimes(indeksy)-samples_skip)*samples_frame);

        if ~isempty(czasy)
            N0=zeros(length(czasy),3);
            N0(:,1) = czasy;
            N0(:,2) = i;
            N0(:,3) = 1;
            N_tmp=[N_tmp' N0']';
        end
    end
end
[t1,ind]=sort(N_tmp(:,1));
N=N_tmp(ind,:);

hold off

xrange = get(graph_axes,'Xlim');
yrange = get(graph_axes,'Ylim');
set(graph_axes,'XLimMode','manual');
set(graph_axes,'YLimMode','manual');
set(graph_axes,'Xlim',xrange);
set(graph_axes,'Ylim',yrange);

if double_mode
    subplot(graph2_axes);
    hold on

    M2=zeros(length(idList2),6);
    N_tmp=[];
    elli2 = zeros(length(idList2),1);
    for i=1:length(idList2)
        neuronID=idList2(i);
        M2(i,1) = i;
        M2(i,2) = paramsFile.getDoubleCell(neuronID(1), 'x0');
        M2(i,3) = paramsFile.getDoubleCell(neuronID(1), 'y0');
        M2(i,4) = paramsFile.getDoubleCell(neuronID(1), 'SigmaX');
        M2(i,5) = paramsFile.getDoubleCell(neuronID(1), 'SigmaY');
        M2(i,6) = paramsFile.getDoubleCell(neuronID(1), 'Theta');

        x1 = M2(i,4)*sin(2*pi*t);
        y1 = M2(i,5)*cos(2*pi*t);
        x2 = x1*cos(M2(i,6)) + y1*sin(M2(i,6)) + M2(i,2);
        y2 = - x1*sin(M2(i,6)) + y1*cos(M2(i,6)) + M2(i,3);
        elli2(i) = fill(y2,x2,'w','LineWidth',1);
      
        if neuronFile.containsID(neuronID)
            spikeTimes = neuronFile.getSpikeTimes(neuronID)';
            indeksy = find(spikeTimes > samples_skip & spikeTimes < (samples_skip + samples_show)); %skip time_offset ms and visualize time_dur ms
            czasy=round((spikeTimes(indeksy)-samples_skip)*samples_frame);

            if ~isempty(czasy)
                N0=zeros(length(czasy),3);
                N0(:,1) = czasy;
                N0(:,2) = i;
                N0(:,3) = 2;
                N_tmp=[N_tmp' N0']';
            end
        end
    end
    [t1,ind]=sort(N_tmp(:,1));
    N2=N_tmp(ind,:);

    hold off

    % Set common axes limits
    xrange = get(graph_axes,'Xlim');
    yrange = get(graph_axes,'Ylim');
    x2range = get(graph2_axes,'Xlim');
    y2range = get(graph2_axes,'Ylim');
    set(graph_axes,'XLimMode','manual');
    set(graph_axes,'YLimMode','manual');
    set(graph2_axes,'XLimMode','manual');
    set(graph2_axes,'YLimMode','manual');
    set(graph_axes,'Xlim',[min(xrange(1),x2range(1)), max(xrange(2),x2range(2))]);
    set(graph_axes,'Ylim',[min(yrange(1),y2range(1)), max(yrange(2),y2range(2))]);
    set(graph2_axes,'Xlim',[min(xrange(1),x2range(1)), max(xrange(2),x2range(2))]);
    set(graph2_axes,'Ylim',[min(yrange(1),y2range(1)), max(yrange(2),y2range(2))]);
end


% Rysowanie dla zerowej klatki
n = 1; %numeruje zdarzenia
light = zeros(length(idList),1); %jak dlugo ma sie swiecic (0 - wylaczony, n > 0 - wylacz po n klatkach)
light = light - 2;
while N(n,1) == 0
    light(N(n,2)) = time_light; %when should light off
    set(elli(N(n,2)),'Visible','on');
    set(elli(N(n,2)),'FaceColor',color(1,:));
    n = n+1;
end
if double_mode
    n2 = 1;
    light2 = zeros(length(idList2),1); %jak dlugo ma sie swiecic (0 - wylaczony, n > 0 - wylacz po n klatkach)
    light2 = light2 - 2;
    while N2(n2,1) == 0
        light2(N2(n2,2)) = time_light; %when should light off
        set(elli2(N2(n2,2)),'Visible','on');
        set(elli2(N2(n2,2)),'FaceColor',color(2,:));
        n2 = n2+1;
    end
end
% Zapis do GIFa
orig_mode = get(fig, 'PaperPositionMode');
set(fig, 'PaperPositionMode', 'auto');
img = hardcopy(fig, '-Dopengl', '-r0');
% Restore figure to original state
set(fig, 'PaperPositionMode', orig_mode); % end
map = [0 0 0; alpha; 1 1 1; 1 0 0; 0 1 0];
A = rgb2ind(img,map,'nodither');
imwrite(A,map,[gif_path '.generating'],'gif','LoopCount',inf,'DelayTime',gif_time_step,'ScreenSize',[gif_width gif_height]);

set(fig,'color',alpha);
set(label,'Visible','off');
set(label2,'Visible','off');

set(graph_axes,'Visible','off');
prop_value = cell(length(elli),1);
prop_value(:) = {'off'};
set(elli,{'Visible'},prop_value);
org_chil = get(graph_axes,'Children');

if double_mode
    set(graph2_axes,'Visible','off');
    prop2_value = cell(length(elli2),1);
    prop2_value(:) = {'off'};
    set(elli2,{'Visible'},prop2_value);
    org_chil2 = get(graph2_axes,'Children');
end

for i = 1:floor(time_dur/time_res)
    subplot(graph_axes);
    hold on
    if n <= length(N) && N(n,1) == i
        while n <= length(N) && N(n,1) == i
            light(N(n,2)) = time_light; %when should light off
            n = n+1;
        end
    end
    
    [~,index] = sort(light(end:-1:1),'descend');
    set(graph_axes,'Children',org_chil(index));
    
    light_off = find(light == 0); %masz je teraz wylaczyc
    if ~isempty(light_off)
        for j = 1:length(light_off)
            set(elli(light_off(j)),'Visible','on');
            set(elli(light_off(j)),'FaceColor','w');
        end
    end
    light_alpha = find(light==-1 | (light==time_light-1 & light > 0)); %masz je teraz znowu zrobic przezroczyste
    if ~isempty(light_alpha)
        for j = 1:length(light_alpha)
            set(elli(light_alpha(j)),'Visible','off');
        end
    end
    light_on = find(light == time_light);
    if ~isempty(light_on)
        for j = 1:length(light_on)
            set(elli(light_on(j)),'Visible','on');
            set(elli(light_on(j)),'FaceColor',color(1,:));
        end
    end
    
    light = light - 1;
    light(light < -2) = -2;
    hold off
    
    if double_mode
        subplot(graph2_axes);
        hold on
        if n2 <= length(N2) && N2(n2,1) == i
            while n2 <= length(N2) && N2(n2,1) == i
                light2(N2(n2,2)) = time_light; %when should light off
                n2 = n2+1;
            end
        end

        [~,index] = sort(light2(end:-1:1),'descend');
        set(graph2_axes,'Children',org_chil2(index));
        light_off = find(light2==0); %masz je teraz wylaczyc (narysuj biale, zeby przykryc kolor)
        if ~isempty(light_off)
            for j = 1:length(light_off)
                set(elli2(light_off(j)),'Visible','on');
                set(elli2(light_off(j)),'FaceColor','w');
            end
        end
        light_alpha = find(light2==-1 | (light2==time_light-1 & light2 > 0)); %masz je teraz znowu zrobic przezroczyste
        if ~isempty(light_alpha)
            for j = 1:length(light_alpha)
                set(elli2(light_alpha(j)),'Visible','off');
            end
        end
        light_on = find(light2 == time_light);
        if ~isempty(light_on)
            for j = 1:length(light_on)
                set(elli2(light_on(j)),'Visible','on');
                set(elli2(light_on(j)),'FaceColor',color(2,:));
            end
        end

        light2 = light2 - 1;
        light2(light2 < -2) = -2;
        hold off
    end
    
    subplot(full_axes);
    set(time_text,'String',[num2str(time_offset_sec + i/400,'%1.3f') ' s']);
    
    img = hardcopy(fig, '-Dopengl', '-r0');
    A = rgb2ind(img,map,'nodither');
    imwrite(A,map,[gif_path '.generating'],'gif','WriteMode','append','DelayTime',gif_time_step,'DisposalMethod','leaveInPlace','TransparentColor',1);
    
    waitbar(i/floor(time_dur/time_res),w);
end
movefile([gif_path '.generating'],[gif_path '.gif'],'f');
close(w);
msg = 'The animation has been successfuly generated!';
msgbox(msg,'Information');
