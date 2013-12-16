function clusters_onefile(params)

    % Saves the params files, 
    % and other statistiscs on the agents level, clusters are not computed
    
    folderName = params.folderName;
    simName = params.simName;
    fileName = params.fileName;
    DUMP = params.DUMP;
    DUMP_RATE = params.DUMP_RATE;
    PLOTS = params.PLOTS;
    CLU_CUTOFF = params.CLU_CUTOFF;
    outDir = params.outDirClusters;
    
    path = [folderName simName fileName];
    load(path);
    
    colnorm = @(X,P) sum(abs(X).^P,1).^(1/P);
  
    % Date and Time
    mytimestamp = datestr ( datevec ( now ), 0 );
    
    v = dump.agentsv;
    pos = dump.agents;
    truth = dump.truth;
    run = dump.run;
    simnameidx = dump.sim;
    
    nIter = size(pos,3);
    nAgents = size(pos, 2);
 
    % GLOBAL variables.
    global_count_sum = zeros(nIter,1);
    global_count_sumsquared = zeros(nIter,1);
    global_maxsize_sum = zeros(nIter,1);
    global_maxsize_sumsquared = zeros(nIter,1);
    global_meansize_sum = zeros(nIter,1);
    global_meansize_sumsquared = zeros(nIter,1); 
    global_speed_sum = zeros(nIter,1);
    global_speed_sumsquared = zeros(nIter,1);
    global_movs_sum = zeros(nIter,1);
    global_movs_sumsquared = zeros(nIter,1);
    global_fromtruth_sum = zeros(nIter,1);
    global_fromtruth_sumsquared = zeros(nIter,1);
    global_bigc_pairwise_pdist_sum = zeros(nIter,1);
    global_bigc_pairwise_pdist_sumsquared = zeros(nIter,1);
    
    if (DUMP)
   
        dataFileName = [outDir 'clusters_macro_' fileName '.csv'];
        fidClustersMacro = fopen(dataFileName, 'w');
        
        dataFileName = [outDir 'clusters_micro_' fileName '.csv'];
        fidClustersMicro = fopen(dataFileName, 'w');     
    end
    
    for i = 1:nIter
       
        % Movement.
        if (i == 1)
            movs = zeros(nAgents, 1);
        else
            movs = colnorm(pos(:,:,i) - pos(:,:,i-1), 2);
        end
        
        
        % If agents are too clustered matlab freezes on LINKAGE
        % calculation. We keep the values from last time.            
        if (pdist(pos(:,:,i)) > 1e-14)
        
            try
                % Z the results of HCLUST 
                % T the cluster of each agent
                % nClusters the number of clusters            
                [Z, T, nClusters] = clusterize(pos(:,:,i));
                
                [d, clusters_size, clusters_fromtruth, clusters_speed, clusters_move] = ...
                cluster_stats(T, nClusters, truth, pos(:,:,i), v(:,:,i), movs);

                max_cluster_size = max(clusters_size);
                bigGroupId = find(clusters_size == max_cluster_size);
                if (length(bigGroupId) > 1) 
                    bigGroupId = bigGroupId(1);
                end

                mean_cluster_size = mean(clusters_size);
                sd_cluster_size = std(clusters_size);

                bigc_agents = pos(:, T == bigGroupId, i);

                pairwise_dist_bigc = pdist(bigc_agents', 'euclidean');
                pairwise_dist_bigc_mean = mean(pairwise_dist_bigc);
                pairwise_dist_bigc_sd = std(pairwise_dist_bigc);

                mean_cluster_fromtruth = mean(clusters_fromtruth);
                sd_cluster_fromtruth = std(clusters_fromtruth);

                mean_cluster_speed = mean(clusters_speed);
                sd_cluster_speed = std(clusters_speed);
                mean_cluster_move = mean(clusters_move);
                sd_cluster_move = std(clusters_move);


            catch err

                err
                sprintf('Error at iter: %i, fileName: %s', i, fileName)

                max_cluster_size = 0;
                mean_cluster_size = 0;
                sd_cluster_size = 0;

                mean_cluster_fromtruth = 0;
                sd_cluster_fromtruth = 0;

                mean_cluster_speed = 0;
                sd_cluster_speed = 0;
                mean_cluster_move = 0;            
                sd_cluster_move = 0;

                % pdist
                pairwise_dist_bigc_mean = 0;
                pairwise_dist_bigc_sd = 0;

                % micro.
                clusters_size = 0;
                clusters_fromtruth = 0;
                clusters_speed = 0;
                clusters_move = 0;

            end


        else                
            sprintf('Too clustered at iter: %i, fileName: %s', i, fileName)
        end
            
        % GLOBAL statistics: will be used by aggregate function.
        global_count_sum(i) = nClusters;
        global_count_sumsquared(i) = nClusters^2;
        global_meansize_sum(i) = mean_cluster_size;
        global_meansize_sumsquared(i) = mean_cluster_size^2;
        global_maxsize_sum(i) = max_cluster_size;
        global_maxsize_sumsquared(i) = max_cluster_size^2;
        global_speed_sum(i) = mean_cluster_speed;
        global_speed_sumsquared(i) = mean_cluster_speed^2;
        global_movs_sum(i) = mean_cluster_move;
        global_movs_sumsquared(i) = mean_cluster_move^2;               
        global_fromtruth_sum(i) = mean_cluster_fromtruth;
        global_fromtruth_sumsquared(i) = mean_cluster_fromtruth^2;
        global_bigc_pairwise_pdist_sum(i) = pairwise_dist_bigc_mean;
        global_bigc_pairwise_pdist_sumsquared(i) = pairwise_dist_bigc_mean^2;

        
        if (DUMP)

             % SAVING ONLY EVERY X ITERATIONS        
             if (mod(i, DUMP_RATE) == 0)

                 stepData = struct(...
                     'simnameidx', simnameidx, ...
                     'run', dump.run, ...
                     'cluster_count', nClusters, ... % count 
                     'mean_cluster_size', mean_cluster_size, ... % size 
                     'sd_cluster_size', sd_cluster_size, ... 
                     'max_cluster_size', mean_cluster_size, ...                      
                     'mean_cluster_speed', mean_cluster_speed, ... % speed
                     'sd_cluster_speed', sd_cluster_speed, ...
                     'mean_cluster_move', mean_cluster_move, ... % move
                     'sd_cluster_move', sd_cluster_move, ...
                     'mean_from_truth', mean_cluster_fromtruth, ... % truth
                     'sd_from_truth', sd_cluster_fromtruth, ...
                     'bigc_pairwise_pdist', pairwise_dist_bigc_mean, ... % bigc dist
                     'bigc_pairwise_pdist_sd', pairwise_dist_bigc_sd, ... % bigc dist
                     'clusters_speed', clusters_speed, ... %micro
                     'clusters_move', clusters_move, ... %micro
                     'clusters_size', clusters_size, ... %micro
                     'clusters_fromtruth', clusters_fromtruth ... %micro
                );

                % 1 Line
                clu_macro_string = csv_format_row_clusters_macro_new(...
                    stepData, simName, i);
                fprintf(fidClustersMacro,'%s\n', clu_macro_string);

                % SAVING ONLY CLUSTERS of SIZE > CUTOFF
                idxs = find(stepData.clusters_size > CLU_CUTOFF);

                % Multiple Lines
                for jj = 1 : length(idxs)
                    clu_micro_string = csv_format_row_clusters_micro( ...
                        stepData, simName, i, idxs(jj));
                    fprintf(fidClustersMicro,'%s\n', clu_micro_string);   
                end
                
             end

        end           
      
        if (PLOTS)
            % something.
        end
        
    end

    save([ outDir 'sums_' fileName ], ...
                            'global_count_sum', ...
                            'global_count_sumsquared', ...
                            'global_meansize_sum', ...
                            'global_meansize_sumsquared', ...
                            'global_maxsize_sum', ...
                            'global_maxsize_sumsquared', ...
                            'global_speed_sum', ...
                            'global_speed_sumsquared', ...
                            'global_movs_sum', ...
                            'global_movs_sumsquared', ...
                            'global_fromtruth_sum', ...
                            'global_fromtruth_sumsquared', ...
                            'global_bigc_pairwise_pdist_sum', ...
                            'global_bigc_pairwise_pdist_sumsquared' ...
    );
 
    if (DUMP)        
        fclose(fidClustersMacro);
        fclose(fidClustersMicro);
    end
end

%%
