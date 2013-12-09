function temporal_analysis(params)

    % Saves the params files, 
    % and other statistiscs on the agents level, clusters are not computed
    
    folderName = params.folderName;
    fileName = params.fileName;
    CSV_DUMP = params.DUMP;
    DUMP_RATE = params.DUMP_RATE;
    PLOTS = params.PLOTS;
    CLU_CUTOFF = params.CLU_CUTOFF;
    outDir = params.outDirClusters;
    
    path = [folderName fileName];
    load(path);
     
    %% Param
             
    headers_clusters_macro = {
        'simname', ...
        'simcount', ...
        'run', ...
        't', ...
        'count', ...
        'size.max', ...
        'size.avg', ...
        'size.sd', ...
        'coverage', ...
        'coverage.cum', ...
        'speed.avg', ...
        'speed.sd', ...
        'move.avg', ...
        'move.sd', ...
        'fromtruth.avg', ...
        'fromtruth.sd' ...
    };

    headers_clusters_micro = {
        'simname', ...
        'simcount', ...
        'run', ...
        't', ...
        'size', ...
        'speed', ...
        'move', ...
        'fromtruth'...
    };
    
    % Date and Time
    mytimestamp = datestr ( datevec ( now ), 0 );
    
    v = dump.agentsv;
    pos = dump.agents;
    truth = dump.truth;
    run = dump.run;
    simName = dump.name;
    simnameidx = dump.sim;

    % GLOBAL variables
    global_count_sum = zeros(nIter,1);
    global_count_sumsquared = zeros(nIter,1);
    global_maxsize_sum = zeros(nIter,1);
    global_maxsize_sumsquared = zeros(nIter,1);
    global_meansize_sum = zeros(nIter,1);
    global_meansize_sumsquared = zeros(nIter,1); 
    global_speed_sum = zeros(nIter,1);
    global_speed_sumsquared = zeros(nIter,1);
    global_move_sum = zeros(nIter,1);
    global_move_sumsquared = zeros(nIter,1);
    global_fromtruth_sum = zeros(nIter,1);
    global_fromtruth_sumsquared = zeros(nIter,1);
    global_bigc_pairwise_pdist_sum = zeros(nIter,1);
    global_bigc_pairwise_pdist_sumsquared = zeros(nIter,1);
    
    if (CSV_DUMP)

        dataFileName = [dumpDir 'clusters_macro.csv'];
        write_csv_headers(dataFileName, headers_clusters_macro);
        fidClustersMacro = fopen(dataFileName,'a');
        
        dataFileName = [dumpDir 'clusters_micro.csv'];
        write_csv_headers(dataFileName, headers_clusters_micro);
        fidClustersMicro = fopen(dataFileName,'a');                
       
        param_string = csv_format_row_params(simName, simnameidx, ...
            run, mytimestamp, dump.parameters, truth);
        
        % append the param string to the file
        fprintf(fidParam,'%s\n', param_string);
    end


    for i = 1:nIter
       
        try
            % Z the results of HCLUST 
            % T the cluster of each agent
            % C the number of clusters            
            [Z, T, C] = clusterize(pos(:,:,i));
            cluster_count = C;

            [d, Gc, AvgGDist, avgGroupSpeed, avgGroupMove] = ...
                cluster_stats(T, truth, pos(:,:,i), v(:,:,i), movs);

            bigGroupId = Gc == max(Gc);
            max_cluster_size = Gc(bigGroupId);
            mean_cluster_size = mean(Gc);
            sd_cluster_size = std(Gc);
        
            bigc_agents = pos(:, T == bigGroupId, i);
            
            pairwise_dist_bigc = pdist(bigc_agents', 'euclidean');
            pairwise_dist_bigc_mean = mean(bigc_agents);
            pairwise_dist_bigc_sd = std(bigc_agents);
                      
            mean_cluster_fromtruth = mean(AvgGDist);
            sd_cluster_fromtruth = std(AvgGDist);

            sd_cluster_speed = std(avgGroupSpeed);
            sd_cluster_move = std(avgGroupMove);

            clusters_size = Gc;
            clusters_fromtruth = AvgGDist;
            clusters_speed = avgGroupSpeed;
            clusters_move = avgGroupMove;

        catch err

            err
            sprintf('Error at iter: %i, fileIdx: %i, fileName: %s', i, validFileIdx, fileName)

            mean_cluster_size = 0;
            sd_cluster_size = 0;

            max_cluster_size = 0;
            
            mean_cluster_fromtruth = 0;
            sd_cluster_fromtruth = 0;

            sd_cluster_speed = 0;
            sd_cluster_move = 0;

            clusters_size = 0;
            clusters_fromtruth = 0;
            clusters_speed = 0;
            clusters_move = 0;

        end

        
        % GLOBAL statistics: will be used by aggregate function.
        global_count_sum(i) = C;
        global_count_sumsquared(i) = C^2;
        global_meansize_sum(i) = mean_cluster_size;
        global_meansize_sumsquared(i) = mean_cluster_size^2;
        global_maxsize_sum(i) = max_cluster_size;
        global_maxsize_sumsquared(i) = max_cluster_size^2;
        global_speed_sum(i) = mean_cluster_speed;
        global_speed_sumsquared(i) = mean_cluster_speed^2;
        global_move_sum(i) = mean_cluster_move;
        global_move_sumsquared(i) = mean_cluster_move^2;               
        global_fromtruth_sum(i) = mean_cluster_fromtruth;
        global_fromtruth_sumsquared(i) = mean_cluster_fromtruth^2;
        global_bigc_pairwise_pdist_sum(i) = pairwise_dist_bigc_mean;
        global_bigc_pairwise_pdist_sumsquared(i) = pairwise_dist_bigc_mean^2;

        
        if (CSV_DUMP)

             % SAVING ONLY EVERY X ITERATIONS        
             if (mod(i, DUMP_RATE) == 0)

                 stepData = struct(...
                     'simnameidx', simnameidx, ...
                     'run', dump.run, ...
                     'cluster_count', cluster_count, ... 
                     'mean_cluster_size', mean_cluster_size, ... 
                     'sd_cluster_size', sd_cluster_size, ... 
                     'max_cluster_size', mean_cluster_size, ...                      
                     'mean_cluster_speed', mean_cluster_speed, ...
                     'mean_cluster_move', mean_cluster_move, ...
                     'sd_cluster_move', sd_cluster_move, ...
                     'sd_cluster_speed', sd_cluster_speed, ...                  
                     'clusters_speed', clusters_speed, ... %micro
                     'clusters_move', clusters_move, ... %micro
                     'clusters_size', clusters_size, ... %micro
                     'clusters_fromtruth', clusters_fromtruth ... %micro
                );

                % 1 Line
                clu_macro_string = csv_format_row_clusters_macro(stepData, simName, i);
                fprintf(fidClustersMacro,'%s\n', clu_macro_string);

                % SAVING ONLY CLUSTERS of SIZE > CUTOFF
                idxs = find(stepData.clusters_size > CLU_CUTOFF);

                % Multiple Lines
                for jj=1:length(idxs)
                    clu_micro_string = csv_format_row_clusters_micro(stepData, simName, i, idxs(jj)); 
                    fprintf(fidClustersMicro,'%s\n', clu_micro_string);   
                end

             end

        end           
      
        if (PLOTS)
            % something.
        end
        
    end

    save([ outDir 'clusters_' filename], ...
                                    'global_count_sum', ...
                                    'global_count_sumsquared', ...
                                    'global_meansize_sum', ...
                                    'global_meansize_sumsquared', ...
                                    'global_maxsize_sum', ...
                                    'global_maxsize_sumsquared', ...
                                    'global_speed_sum', ...
                                    'global_speed_sumsquared', ...
                                    'global_move_sum', ...
                                    'global_move_sumsquared', ...
                                    'global_fromtruth_sum', ...
                                    'global_fromtruth_sumsquared' ...
    );
    
 
    
    
    
end

%%




