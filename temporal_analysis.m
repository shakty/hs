function temporal_analysis( DUMPDIR, simName, PRECISION, CSV_DUMP, PLOTS)

    
    %% Param
             
    headers_clusters_macro = {
        'simname', ...
        'simcount', ...
        'run', ...
        't', ...
        'count', ...
        'coverage', ...
        'coverage.cum', ...
        'speed.avg', ...
        'speed.sd', ...
        'move.avg', ...
        'move.sd', ...
        'size.avg', ...
        'size.sd', ...
        'fromtruth.avg', ...
        'fromtruth.sd'};
    
    headers_clusters_micro = {
        'simname', ...
        'simcount', ...
        'run', ...
        't', ...
        'size', ...
        'speed', ...
        'move', ...
        'fromtruth'};

    headers_params = {
        'simname', ...
        'simcount', ...
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
    
    
    TOTAL_CELLS = PRECISION^2;

    
    colnorm = @(X,P) sum(abs(X).^P,1).^(1/P);

    dumpDir = [DUMPDIR simName '/'];
    
    files = dir(dumpDir);
    fileIndex = find(~[files.isdir]);

    if (isempty(fileIndex))
        error('Invalid Directory Selected');
    end

    % Date and Time
    mytimestamp = datestr ( datevec ( now ), 0 );

    % this can be longer than needed, but does not matter
    simnameidxs = zeros(length(fileIndex),1);
    runs = zeros(length(fileIndex),1);
    
    if (CSV_DUMP)

        % params (for both)
        paramFileName = [dumpDir 'params.csv'];
        write_csv_headers(paramFileName, headers_params);
        fidParam = fopen(paramFileName,'a');

        
        dataFileName = [dumpDir 'clusters_macro.csv'];
        write_csv_headers(dataFileName, headers_clusters_macro);
        fidClustersMacro = fopen(dataFileName,'a');
        
        dataFileName = [dumpDir 'clusters_micro.csv'];
        write_csv_headers(dataFileName, headers_clusters_micro);
        fidClustersMicro = fopen(dataFileName,'a');

    end
    

    % Load all parameters matrices in one
    for f = 1:length(fileIndex)

        append = files(fileIndex(f)).name;
        fileName = [dumpDir, append];

        % We load only .mat
        [PATH,NAME,EXT] = fileparts(fileName);
        if (~strcmpi(EXT,'.mat')) 
            continue;
        end

        simnameidx = strfind(NAME, '-');
        simnameidx = NAME(1:simnameidx-1);
        simnameidxs(f) = simnameidx;
        
        load(fileName);
        
        runs(f) = dump.run;
        
        %dump.parameters
        
        if (CSV_DUMP)
            param_string = csv_format_row_params(simName, simnameidx, dump.run, mytimestamp, dump.parameters, dump.truth);
            % append the param string to the file
            fprintf(fidParam,'%s\n', param_string);
        end
    
        

        v = dump.agentsv;
        pos = dump.agents;
        
        nIter = size(v,3);
        
        % average velocity of agents at time t
        mean_cluster_speed = zeros(1,nIter);
        % Std velocity of clusters at time t
        sd_cluster_speed = zeros(1,nIter);

        % vector of movements of agents at time t
        movs = zeros(size(pos,2),1);
        % average movement from position at time t-1
        mean_cluster_move = zeros(1,nIter);
        % Std movement of clusters at time t-1
        sd_cluster_move = zeros(1,nIter);

        
        % average share of space occupied by agents at time t
        avgcoverage = zeros(1,nIter);
        % cumulative share of space explored by agents at time t
        cumcoverage = zeros(1,nIter);
        % cumulative occupation of each cell of the grid at time t
        cum_coverage_matrix = zeros(PRECISION, PRECISION, 3);
        % cumulative occupation of each cell of the grid at time t
        coverage_matrix = zeros(PRECISION, PRECISION, 3);

        % Counts the number of clusters
        cluster_count = zeros(1,nIter);              
        % Average cluster size
        mean_cluster_size = zeros(1,nIter);          
        % Standard deviations of the size of clusters
        sd_cluster_size = zeros(1,nIter);          	
        % Average distance from truth
        mean_from_truth = zeros(1, nIter);    
        % Standard deviation of distance from truth
        sd_from_truth = zeros(1, nIter);    


        % Following cell arrays containing vector of variable length at each iteration

        % Size of the clusters at time t
        clusters_size = cell(nIter,1);
        % Speed of clusters at time t
        clusters_speed = cell(nIter,1);
        % Spatial displacement of the cluster compared with time t-1
        % It is the average displacement of the agents within in
        clusters_move = cell(nIter,1);
        % Distance from truth at time t
        clusters_fromtruth = cell(nIter,1);

        for i = 1:nIter
            %avg speed
            mean_cluster_speed(i) = mean(colnorm(v(:,:,i),2));
            % avg movement
            if (i>1)
                movs = abs(pos(:,:,i)-pos(:,:,i-1));
                mean_cluster_move(i) = mean(mean(movs));
            end

            %avg share of space
            coverage_matrix(:,:,i) = countAgents(pos(:,:,i), PRECISION);
            avgcoverage(i) = nnz(coverage_matrix(:,:,i)) / TOTAL_CELLS;

            %cum share of space
            if (i==1)
                cum_coverage_matrix(:,:,i) = coverage_matrix(:,:,i);
            else
                cum_coverage_matrix(:,:,i) = coverage_matrix(:,:,i) + cum_coverage_matrix(:,:,i-1);
            end
            cumcoverage(i) = nnz(cum_coverage_matrix(:,:,i)) / TOTAL_CELLS;

            % Z the results of HCLUST 
            % T the cluster of each agent
            % C the number of clusters
            [Z, T, C] = clusterize(pos(:,:,i));
            cluster_count(i) = C;

            [d, Gc, AvgGDist, avgGroupSpeed, avgGroupMove] = cluster_stats(T, dump.truth, pos(:,:,i), v(:,:,i), movs);

            mean_cluster_size(i) = mean(Gc);
            sd_cluster_size(i) = std(Gc);
            mean_from_truth(i) = mean(AvgGDist);
            sd_from_truth(i) = std(AvgGDist);
            sd_cluster_speed(i) = std(avgGroupSpeed);
            sd_cluster_move(i) = std(avgGroupMove);

            clusters_size{i} = Gc;
            clusters_fromtruth{i} = AvgGDist;
            clusters_speed{i} = avgGroupSpeed;
            clusters_move{i} = avgGroupMove;

        end

    end
    
    hold on
    %plot(1:nIter, mean_cluster_speed, 'r')
    %plot(1:nIter, mean_cluster_move, 'b')
    plot(1:nIter, avgcoverage, 'g')
    plot(1:nIter, cumcoverage, 'g')
    hold off



    if (CSV_DUMP)
        for z = 1:nIter                
            cluObj = struct(...
                 'name', simName, ...
                 'simnameidx', simnameidx(z), ...
                 'run', runs(z), ...
                 't', z, ...
                 'mean_cluster_speed', mean_cluster_speed(z), ...
                 'mean_cluster_move', mean_cluster_move(z), ...
                 'sd_cluster_move', sd_cluster_move(z), ...
                 'sd_cluster_speed', sd_cluster_speed(z), ...
                 'cluster_count', cluster_count(z), ...
                 'avgcoverage', avgcoverage(z), ...
                 'cumcoverage', cumcoverage(z), ...
                 'mean_cluster_size', mean_cluster_size(z), ... 
                 'sd_cluster_size', sd_cluster_size(z), ... 
                 'mean_from_truth', mean_from_truth(z), ... 
                 'sd_from_truth', sd_from_truth(z), ...
                 'clusters_sizes', clusters_size{z}, ...
                 'clusters_from_fromtruth', clusters_fromtruth{z}, ...
                 'clusters_speed', clusters_speed{z}, ...
                 'clusters_move', clusters_move{z} ...
            );

            clu_macro_string = csv_format_row_clusters_macro(cluObj);
            fprintf(fidClustersMacro,'%s\n', clu_macro_string);
            
            
            for k = 1:length(clusters_size(z))
                clu_micro_string = csv_format_row_clusters_macro(cluObj); 
                fprintf(fidClustersMicro,'%s\n', clu_micro_string);
            end
            
            
        end
    end
    
    if (CSV_DUMP)
    
        fclose(fidParam);
        fclose(fidCorr);
        fclose(fidClusters);
    end
end

%%




