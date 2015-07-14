%% Start analysis; 
% Script that looks at the activation probability of neurons after
% stimulation along the axon. The purpose here is to study the change in
% activation probability using 2-electrode, opposite polarity stimulation
% with respect to 1-electrode stimulation. 

neuronIds = {'602'; '1171'; '1699' ; '5267'};
[xc,yc] = getElectrodeCoords512(); 
% Load red teal diverging colormap
codebase_path = matlab_code_path; 
cmap = load([codebase_path 'code/projects/electrical_stim/'...
    'resources/redtealcmap.mat']); 

pairs = []; 
ii = 1; 
for n = 1:length(neuronIds)
    neuronId = str2double(neuronIds{n});
    switch neuronId
        case 602
            fpath = '/Volumes/Analysis/2015-04-14-0/data001/';
            fpath2elec = '/Volumes/Analysis/2015-04-14-0/data002/';
        case 1171
            fpath = '/Volumes/Analysis/2015-04-14-0/data001/';
            fpath2elec = '/Volumes/Analysis/2015-04-14-0/data002/';
        case 1699
            fpath = '/Volumes/Analysis/2015-04-14-0/data001/';
            fpath2elec = '/Volumes/Analysis/2015-04-14-0/data002/';
        case 5267
            fpath = '/Volumes/Analysis/2015-04-14-0/data001/';
            fpath2elec = '/Volumes/Analysis/2015-04-14-0/data002/';
        otherwise
            return;
    end
    [eiAmps, allVals, actProb, actAmp] = genActThreshSpatialMaps(fpath,...
        neuronId, 'printElecs',1,'suppressPlots',1);
     
    % actProb = the probability of cell activation at actAmp
    % actAmp, use this amplitude to get the probability with a given 2 elec
    % pattern
    
    % Plot the EI contour with a linear fit + electrodes
    [XI,YI] =  eiContour_wLinFit(eiAmps,'linFitThresh',5);
    title(sprintf('neuron %d',neuronId));
    
    % Show stimulating electrodes
    activatingElecs = find(allVals);
    hold on; scatter(xc(activatingElecs),yc(activatingElecs),350, ...
        allVals(activatingElecs),'filled');
    colormap(cmap.m); caxis([0.5 3.5]); colorbar; 
    
    patternPath = [fpath2elec 'pattern_files/'];
    elecsUsed = getElecsFromPatternFiles(patternPath); 
    [~, biElecThresholds] = genActThreshSpatialMaps(fpath2elec,...
        neuronId, 'assignValToPatternNo',1,'suppressPlots',1);
    for e = 1:size(activatingElecs,2)
        [biElecPatterns,~] = find(activatingElecs(e) == elecsUsed); 
        % Find patterns where primary electrode is negative and secondary
        % is positive
        
        for ee = 1:length(biElecPatterns)
            movieNos = findMovieNos(fpath2elec,biElecPatterns(ee));
            [amps, elecs] = getStimAmps(fpath2elec, ...
                biElecPatterns(ee), movieNos(end));
            primaryElecIdx = find(elecs==activatingElecs(e));
            if amps(primaryElecIdx) < 0 % Only analyze pairs with negative current amplitude on the primary electrode.
%                 biElecThresholds(biElecPatterns(ee))
                if primaryElecIdx == 2
                    elecs = flip(elecs,2); 
                end
                primElecProb = actProb(elecs(1)); %gives the activation probability at the primary electrode
                ampOfInterest = actAmp(elecs(1)); %gives the current at that electrode where the activation happened. 
                % What is the probability of activation of current pattern
                % biElecPatterns(ee) at this stimulation current amplitude?
                try
                    respProb = getRespProb(fpath2elec,biElecPatterns(ee),neuronId,ampOfInterest);
                catch
                    try
                        respProb = getRespProb(fpath2elec,biElecPatterns(ee),1756,ampOfInterest);
                    catch
                        respProb = 0;
                        disp(['warning - did not test pattern ' ...
                            num2str(biElecPatterns(ee)) ' for n' num2str(neuronId)]);
                    end
                end
                 
                % Analyze the current pattern!
                xvals = xc(elecs);
                yvals = yc(elecs);
                hold on; line(xvals,yvals,'LineWidth',2,'Color','k'); %colors(ee,:));
                
                % Find corresponding line segment
                idx1 = find(xvals(1)< YI,1,'first');
                idx2 = find(xvals(1)> YI,1,'last');
                if ~isempty(idx1) && ~isempty(idx2) % dont analyze the array edges
                    xc_line = YI([idx1 idx2]);
                    yc_line = XI([idx1 idx2]);
                    hold on; line(xc_line,yc_line,'LineWidth',2,'Color','k'); %colors(ee,:));
                    
                    a = yc_line(2) - yc_line(1);
                    b = xc_line(2) - xc_line(1);
                    c = xc_line(2)*yc_line(1) - xc_line(1)*yc_line(2);
                    e1_d = abs(a*xvals(1) - b*yvals(1) + c)/sqrt(a^2 + b^2);
                    e2_d = abs(a*xvals(2) - b*yvals(2) + c)/sqrt(a^2 + b^2);
                    d = sqrt(diff(xvals)^2 + diff(yvals)^2);
                    pairs = [pairs; primElecProb respProb];
                    if primElecProb < respProb
                        fprintf(['check EI fit and activation curves '...
                            'for patterns %0.0f, %0.0f neuron: %0.0f\n'],...
                            biElecPatterns(ee),elecs(1),neuronId);
                    end
                    ratio(ii) = e1_d/d;
                    text(xvals(2)+10,yvals(2),sprintf('%0.0f to %0.0f: %0.2f',e1_d,...
                        e2_d,ratio(ii)),'FontSize',16);
%                     fprintf('%0.0f to %0.0f: %0.2f\n',e1_d,e2_d,ratio(ii))
                    ii = ii+1; 
                end
                    
            end
        end
        
        % Find patterns where the secondary electrode crosses the axon. Do
        % this by checking the distances to the EI fit and then finding the
        % minimum. 
    end % End loop through the activating electrodes searching for bielec stim patterns
    
end
colors = lines(1); 
figure; plot(ratio,pairs(:,1) - pairs(:,2),'o','MarkerSize',12,'MarkerFaceColor',colors(1,:)); 
% hold on; line([1 1],[0 1],'LineStyle','--'); 
xlim([0 0.5]); 
xlabel('distance of the primary electrode to the axon / electrode spacing'); 
ylabel('change in activation probability with second electrode'); 

figure; plot([0 -1],pairs','x-','LineWidth',2);
xlabel('relative current through electrode 2'); 
ylabel('activation probability');
title(sprintf('Axonal activation at varied locations\n1- and 2-electrode stimulation at fixed amplitudes'));
