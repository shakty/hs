function [ output_args ] = tocsv( simName, CSV_CLU, CSV_POS, PLOT_POS, PLOT_CLU  )

    %% Param

    headers_pos = {
        'sim', ...
        'run', ...
        't', ...
        'id', ...
        'x', ...
        'y'};

    headers_clusters = {
        'sim', ...
        'run', ...
        't', ...
        'count', ...
        'size.avg', ...
        'size.sd', ...
        'fromtruth.avg', ...
        'fromtruth.sd'};

    headers_params = {
        'sim', ...
        'run', ...
        'timestamp', ...
        't.end', ...
        'dt', ...
        'nagents', ...
        'spacesize',...
        'spacedim', ...
        'alpha', ... % Own velocity
        'R', ...
        'k', ... % Exponent forces
        'A', ... % Attractive force
        'd0', ... 
        'B', ... % Repulsive force
        'd1', ... 
        'tau', ... % Truth strength
        'sigma', ... % Std noise
        'init.vscaling', ...
        'init.nclusters', ...
        'init.clusterradio', ...
        'truth.x', ...
        'truth.y', ...
        'conv'};


    DUMPDIR = 'dump/';


    %% Load Data


    %simName = 'few_big_groups-DIM-vs-ALPHA';

    dumpDir = [DUMPDIR simName '/'];


    files = dir(dumpDir);
    fileIndex = find(~[files.isdir]);

    if (isempty(fileIndex))
        error('Invalid Directory Selected');
    end

    % Date and Time
    mytimestamp = datestr ( datevec ( now ), 0 );


    if (CSV_CLU || CSV_POS)

        % params (for both)
        paramFileName = [dumpDir 'params.csv'];
        write_csv_headers(paramFileName, headers_params);
        fidParam = fopen(paramFileName,'a');

        if (CSV_POS)
            dataFileName = [dumpDir 'positions.csv'];
            write_csv_headers(dataFileName, headers_pos);
            fidData = fopen(dataFileName,'a');
        end

        if (CSV_CLU)
            clusterFileName = [dumpDir 'clusters.csv'];
            write_csv_headers(clusterFileName, headers_clusters);
            fidClusters = fopen(clusterFileName,'a');
        end

    end


    % Load all parameters matrices in one
    for i = 1:length(fileIndex)

        append = files(fileIndex(i)).name;
        fileName = [dumpDir, append];

        % We load only .mat
        [PATH,NAME,EXT] = fileparts(fileName);
        if (~strcmpi(EXT,'.mat')) 
            continue;
        end

        load(fileName);

        roundsIdx = [1:(1 /dump.parameters.dt):length(dump.agents)];
        rounds = dump.agents(:,:,roundsIdx);

        nRounds = size(rounds,3);

        if (CSV_CLU || CSV_POS)
            param_string = csv_format_row_params(simName, dump.run, mytimestamp, dump.parameters, dump.truth, dump.conv);
            % append the param string to the file
            fprintf(fidParam,'%s\n', param_string);
        end


        Ccount = zeros(1,nRounds);              % Counts the number of clusters
        Csize_mean = zeros(1,nRounds);          % Average cluster size
        Csize_sd = zeros(1,nRounds);          	% Standard deviations of the size of clusters

        Cfromtruth_mean = zeros(1, nRounds);    % Average distance from truth
        Cfromtruth_sd = zeros(1, nRounds);      % Standard deviation of distance from truth

        % Vectors variable length
        % Csize = zeros(1, nRounds);            % Vector of the size of clusters
        % Cfromtruth = zeros(1, nRounds);       % 

        for z=1:nRounds

           agentpos = rounds(:,:,z);

           if (PLOT_POS)
                hold on;
                plot(agentpos(1,:),agentpos(2,:),'rx');
                plot(dump.truth(1),dump.truth(2),'go');
                xlim([0,1])
                ylim([0,1])
                hold off;
                pause(0.1)
           end

           [Z, T, C] = clusterize(agentpos);
           Ccount(1,z) = C;

           [d, Gc, AvgGDist] = each_cluster_from_truth(T, dump.truth, agentpos);

           Csize_mean(1,z) = mean(Gc);
           Csize_sd(1,z) = std(Gc);
           Cfromtruth_mean(1,z) = mean(AvgGDist);
           Cfromtruth_sd(1,z) = std(AvgGDist);


           if (CSV_POS)
                for id = 1:size(agentpos,2)         
                    clu_string = csv_format_row_clusters(simName, dump.run, z, id, agentpos(:,id));
                    fprintf(fidClusters,'%s\n', clu_string);
                end
           end

           if (CSV_CLU)

                clu_string = csv_format_row_clusters(simName, dump.run, z, ...
                    C, Csize_mean(1,z), Csize_sd(1,z), Cfromtruth_mean(1,z), ...
                    Cfromtruth_sd(1,z));
                fprintf(fidClusters,'%s\n', clu_string);
           end

        end


        if (PLOT_CLU)
            subplot(2,3,1);
            plot([1:nRounds], Ccount);
            subplot(2,3,2);
            plot([1:nRounds], Csize_mean);
            subplot(2,3,3);
            plot([1:nRounds], Csize_sd);
            subplot(2,3,4);
            plot([1:nRounds], Cfromtruth_mean);
            subplot(2,3,5);
            plot([1:nRounds], Cfromtruth_sd);
            pause(0.5)
        end

    end


    if (CSV_CLU || CSV_POS)

        fclose(fidParam);

        if (CSV_POS)
            fclose(fidData);
        end

        if (CSV_CLU)
             fclose(fidClusters);
        end

    end

end

