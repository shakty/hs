function [ output_args ] = tocsv( simName, agentpos)

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



    % Date and Time
    mytimestamp = datestr ( datevec ( now ), 0 );


   
        % params (for both)
        paramFileName = [dumpDir 'params.csv'];
        write_csv_headers(paramFileName, headers_params);
        fidParam = fopen(paramFileName,'a');

     
        clusterFileName = [dumpDir 'clusters.csv'];
        write_csv_headers(clusterFileName, headers_clusters);
        fidClusters = fopen(clusterFileName,'a');
        




            param_string = csv_format_row_params(simName, dump.run, mytimestamp, dump.parameters, dump.truth, dump.conv);
            % append the param string to the file
            fprintf(fidParam,'%s\n', param_string);
       


        Ccount = zeros(1,nRounds);              % Counts the number of clusters
        Csize_mean = zeros(1,nRounds);          % Average cluster size
        Csize_sd = zeros(1,nRounds);          	% Standard deviations of the size of clusters

        Cfromtruth_mean = zeros(1, nRounds);    % Average distance from truth
        Cfromtruth_sd = zeros(1, nRounds);      % Standard deviation of distance from truth

        % Vectors variable length
        % Csize = zeros(1, nRounds);            % Vector of the size of clusters
        % Cfromtruth = zeros(1, nRounds);       % 

      
           [Z, T, C] = clusterize(agentpos);
           Ccount(1,z) = C;

           [d, Gc, AvgGDist] = each_cluster_from_truth(T, dump.truth, agentpos);

           Csize_mean(1,z) = mean(Gc);
           Csize_sd(1,z) = std(Gc);
           Cfromtruth_mean(1,z) = mean(AvgGDist);
           Cfromtruth_sd(1,z) = std(AvgGDist);


      
        
                clu_string = csv_format_row_clusters(simName, dump.run, z, ...
                    C, Csize_mean(1,z), Csize_sd(1,z), Cfromtruth_mean(1,z), ...
                    Cfromtruth_sd(1,z));
                fprintf(fidClusters,'%s\n', clu_string);
          



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

