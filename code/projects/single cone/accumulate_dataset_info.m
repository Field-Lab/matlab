
% get eccentricity, l/m ratio for reach piece

if ~exist('dinfo','var')
    %nicknames = {'blueberry','kiwi','plantain','peach','grapes','apple','apricot'};
    nicknames = {'apricot'};

    dinfo = struct;

    for nn=1:length(nicknames)
        nickname = nicknames{nn};
        load([single_cone_path 'saved/' nickname '.mat'])

        datarun = load_index(datarun);

        % note nickname
        dinfo(nn).nickname = nickname;

        % get eccentricity info
        dinfo(nn).eccentricity_radius = datarun.piece.eccentricity_radius;
        dinfo(nn).eccentricity_clock = datarun.piece.eccentricity_clock;
        dinfo(nn).eye = datarun.piece.eye;

        % get L/M ratio
        datarun = import_single_cone_data(datarun,nickname);
        dinfo(nn).L_count  = sum(datarun.cones.types == 'L');
        dinfo(nn).M_count  = sum(datarun.cones.types == 'M');
        dinfo(nn).S_count  = sum(datarun.cones.types == 'S');

    end
end

fprintf('\n\n\t%%L\t%%M\t%%S\tL/(L+M)\tt.e.e.\n')
for nn=1:length(nicknames)
    di = dinfo(nn);
    tot_cones =  di.L_count  + di.M_count + di.S_count;
    fprintf('%s\t%4.1f\t%4.1f\t%4.1f\t%4.2f\t',...
        nicknames{nn}(1:4),100*di.L_count/tot_cones,100*di.M_count/tot_cones,100*di.S_count/tot_cones,...
        di.L_count/(di.M_count+di.L_count))
    if di.eccentricity_radius == -1
        fprintf('????\n')
    else
        fprintf('%4.1f\n',di.eccentricity_radius)
    end
end

