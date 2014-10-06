function results = compute_ei_axon_match(axon,ei,positions, varargin)
% compute_ei_axon_match     summarize how well an EI and an axon match
%
% usage:  result = compute_ei_axon_match(axon,ei,positions, varargin)
%
% arguments:     axon - Nx2 matrix of axon path, starting at soma
%                  ei - ExT matrix of EI intensity on each electorde at each time point
%           positions - positions of electrodes
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     result - result of computation
%
%
% optional params, their default values, and what they specify:
%
% method          	'basic'       	how to compare
%                                       'basic' - ???
%
% corrfunc          '@xcorr'        Function for calculating the correlation
%                                   between EIs.  Defaults to local simple
%                                   xcorr function, but can be any other
%                                   2-argument function handle, e.g. @corr
%                                   for linear correlation coefficient or 
%                                   @(x,y) corr(x,y,'type','Spearman') for
%                                   rank correlation.
%
% 2010-01  gauthier
% 2013-11  phli, tweaking, expanding
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('method', 'basic');
p.addParamValue('corrfunc', []);
p.addParamValue('plot', false);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

if isempty(ei)
    error('ei is empty')
end

if isempty(params.corrfunc)
    params.corrfunc = @xcorr;
end


% initialize output struct
results = struct;

switch params.method
    case 'basic'
        % initialize variable noting the relevant electrodes
        keep = true(size(ei,1),1);
        
        % convert EI to just amplitudes
        [eia, peaktimes] = max(abs(ei),[],2);
        
        % Cut soma from EI
%         [~, maxelec] = max(eia);
%         cuttime = peaktimes(maxelec)+2;
%         keep = keep & peaktimes > cuttime;
        
        % compute predicted ei
        [predicted_axon,predicted_soma,axon_dists,soma_dists] = generate_predicted_ei(axon(1:end,:),positions);
        
        switch 3
            case 1 % fit axon and soma
                
                % fit axon and soma
                fit_weights = fminsearch(@(x)sum((  eia - (abs(x(1))*predicted_axon + abs(x(2))*predicted_soma)   ).^2),[1 1]);
                
                % use fit values
                predicted_ei = abs(fit_weights(1))*predicted_axon + abs(fit_weights(2))*predicted_soma;
                
                % only keep points where axon > soma
                keep = keep & (abs(fit_weights(1))*predicted_axon > abs(fit_weights(2))*predicted_soma )';
                
            case 2 % ignore some electrodes
                predicted_ei = predicted_axon;


                % exclude axons near soma
                %keep = keep & (soma_dists' > 100);
                
                % exclude electrodes with largest amplitude (probably somatic spikes)
                %[junk,largest_amplitude] = sort(eia,'descend');
                %keep(largest_amplitude(1:3))=false;

                % exclude distant electrodes
                %keep=keep & (axon_dists<300);

                % exclude electodes with small signal
                %keep = keep & (eia > 5)';
                %thresh = 5*robust_std(eia) + median(eia);
                %keep = keep & (eia > thresh)';

                %keep = keep & (positions(:,1) > -500)';

            case 3
                predicted_ei = predicted_axon + 20*predicted_soma;
        end
        
        
        %fprintf('using %d of %d electrodes\n',sum(keep),length(keep))
        
        
        % compute correlation
        if any(keep)
            results.corr = params.corrfunc(predicted_ei(keep), eia(keep));
        else
            results.corr = NaN;
        end

        if params.plot % plot
            figure(params.plot);clf;
            %axis image;axis xy;colormap gray;hold on;

            subplot(211)
            plot(axon(:,1),axon(:,2),'Color','r');hold on
            % ei in blue
            plot_ei_(eia(keep), positions(keep,:), 0,'alpha',0,'pretty_axes',0,'max_scale',45,'pos_color','b','scale',1,'cutoff',-1)
            % unused electrodes in black
            %plot_ei_(~eia(keep), positions(~keep,:), 0,'alpha',0,'pretty_axes',0,'max_scale',45,'pos_color','k','scale',1,'cutoff',-1)
            
            %plot_ei_(predicted_ei(keep), positions(keep,:), 0,'alpha',0,'pretty_axes',0,'max_scale',45,'pos_color','g','scale',1,'cutoff',-1)
            plot_ei_(predicted_axon(keep), positions(keep,:), 0,'alpha',0,'pretty_axes',0,'max_scale',45,'pos_color','r','scale',1,'cutoff',-1)
            plot_ei_(predicted_soma(keep), positions(keep,:), 0,'alpha',0,'pretty_axes',0,'max_scale',45,'pos_color','m','scale',1,'cutoff',-1)

            set(gca,'ydir','reverse','xdir','reverse')
            
            %figure(31);clf;
            subplot(223)
            plot(axon_dists(keep),eia(keep),'.');xlabel('distance from axon');ylabel('ei amplitude')
            subplot(224)
            plot(eia(keep),predicted_ei(keep),'.');xlabel('ei');ylabel('fake ei')
            title(num2str(results.corr))
        end
end

function xc = xcorr(pred_ei, ei)
xc = dot(pred_ei, ei) / norm(pred_ei) / norm(ei);