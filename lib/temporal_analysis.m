function temporal_analysis( DUMPDIR, simName, PRECISION, CLU_CUTOFF, CSV_DUMP, DUMP_RATE, PLOTS)

    
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
        'fromtruth.sd' ...
        };
    
    
    headers_clusters_macro_avg = {
        'simname', ...
        'simcount', ...
        'run', ...
        't', ...
        'count.avg', ...
        'count.sd', ...
        'count.se', ...
        'count.ci', ...
        'coverage.avg', ...
        'coverage.sd', ...
        'coverage.se', ...
        'coverage.ci', ...
        'coverage.cum.avg', ...
        'coverage.cum.sd', ...
        'coverage.cum.se', ...
        'coverage.cum.ci', ...
        'speed.avg', ...
        'speed.sd', ...
        'speed.se', ...
        'speed.ci', ...
        'move.avg', ...
        'move.sd', ...
        'move.se', ...
        'move.ci', ...
        'size.avg', ...
        'size.sd', ...
        'size.se', ...
        'size.ci', ...
        'fromtruth.avg', ...
        'fromtruth.sd', ...
        'fromtruth.se', ...
        'fromtruth.ci' ...
        };
    
        
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
        'noisetype', ...
        'attrtype' ...
        };
    
    
    % CI
    CI_INT = 0.95/2 + 0.5;

    
    % TOTAL CELLS
    TOTAL_CELLS = PRECISION^2;

    
    colnorm = @(X,P) sum(abs(X).^P,1).^(1/P);

    dumpDir = [DUMPDIR simName '/'];
    
    files = dir(dumpDir);
    fileIndex = find(~[files.isdir]);

    if (isempty(fileIndex))
        error('Invalid Directory Selected');
    end

    validFileIdx = 0;
    summaryObj = {};
    
    
    % Retrieving global values for the set of simulations from the first one
    
    fileName = [dumpDir '1-1.mat'];    
    load(fileName);
    v = dump.agentsv;
    % the total number of time steps per run (assumed constant) and dump rate
    nIter = size(v,3);
    idxsIters = find(mod(1:nIter, DUMP_RATE) == 0);
    % number of files
    nFiles = length(fileIndex);
    
    % GLOBAL variables
    global_count_sum = zeros(nIter,1);
    global_count_sumsquared = zeros(nIter,1);
    global_coverage_sum = zeros(nIter,1);
    global_coverage_sumsquared = zeros(nIter,1);
    global_coverage_cum_sum = zeros(nIter,1);
    global_coverage_cum_sumsquared = zeros(nIter,1);
    global_speed_sum = zeros(nIter,1);
    global_speed_sumsquared = zeros(nIter,1);
    global_move_sum = zeros(nIter,1);
    global_move_sumsquared = zeros(nIter,1);
    global_size_sum = zeros(nIter,1);
    global_size_sumsquared = zeros(nIter,1);
    global_fromtruth_sum = zeros(nIter,1);
    global_fromtruth_sumsquared = zeros(nIter,1);
    
    % Date and Time
    mytimestamp = datestr ( datevec ( now ), 0 );
    
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
                
        dataFileName = [dumpDir 'clusters_macro_avg.csv'];
        write_csv_headers(dataFileName, headers_clusters_macro_avg);
        fidClustersMacroAvg = fopen(dataFileName,'a');
       
    end
    
    % Load all parameters matrices in one
    for f = 1:nFiles

        append = files(fileIndex(f)).name;
        fileName = [dumpDir, append];

        % We load only .mat
        [PATH,NAME,EXT] = fileparts(fileName);
        if (~strcmpi(EXT,'.mat') ...
            || strcmp(NAME, 'temporalysis') == 1 ...
            || strcmp(NAME, 'global_count_sum') == 1 ...
            || strcmp(NAME, 'global_count_sumsquared') == 1 )
            continue;
        end

        simnameidx = strfind(NAME, '-');
        simnameidx = str2double(NAME(1:simnameidx-1));
        
        
        load(fileName);
        validFileIdx = validFileIdx + 1;
        
        %dump.parameters
        
        if (CSV_DUMP)
            param_string = csv_format_row_params(simName, simnameidx, dump.run, mytimestamp, dump.parameters, dump.truth);
            % append the param string to the file
            fprintf(fidParam,'%s\n', param_string);
        end
    
        v = dump.agentsv;
        pos = dump.agents;
        
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
        mean_cluster_fromtruth = zeros(1, nIter);    
        % Standard deviation of distance from truth
        sd_cluster_fromtruth = zeros(1, nIter);    


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
            
            mean_cluster_fromtruth(i) = mean(AvgGDist);
            sd_cluster_fromtruth(i) = std(AvgGDist);
            
            sd_cluster_speed(i) = std(avgGroupSpeed);
            sd_cluster_move(i) = std(avgGroupMove);

            clusters_size{i} = Gc;
            clusters_fromtruth{i} = AvgGDist;
            clusters_speed{i} = avgGroupSpeed;
            clusters_move{i} = avgGroupMove;
    
                
            % SUMMING UP AVG statistics
            if (i == 3)
                test_c(validFileIdx) = C;
            end
            global_count_sum(i) = global_count_sum(i) + C;
            global_count_sumsquared(i) = global_count_sumsquared(i) + C^2;
            global_coverage_sum(i) = global_coverage_sum(i) + avgcoverage(i);
            global_coverage_sumsquared(i) = global_coverage_sumsquared(i) + avgcoverage(i)^2;
            global_coverage_cum_sum(i) = global_coverage_cum_sum(i) + cumcoverage(i);
            global_coverage_cum_sumsquared(i) = global_coverage_cum_sumsquared(i) + cumcoverage(i)^2;
            global_speed_sum(i) = global_speed_sum(i) + mean_cluster_speed(i);
            global_speed_sumsquared(i) = global_speed_sumsquared(i) + mean_cluster_speed(i)^2;
            global_move_sum(i) = global_move_sum(i) + mean_cluster_move(i);
            global_move_sumsquared(i) =  global_move_sumsquared(i) + mean_cluster_move(i)^2;
            global_size_sum(i) = global_size_sum(i) + mean_cluster_size(i);
            global_size_sumsquared(i) = global_size_sumsquared(i) + mean_cluster_size(i)^2;
            global_fromtruth_sum(i) = global_fromtruth_sum(i) +  mean_cluster_fromtruth(i);
            global_fromtruth_sumsquared(i) = global_fromtruth_sumsquared(i) +  mean_cluster_fromtruth(i)^2;

        end
        
        
        summaryObj{validFileIdx} = struct(...
             'simnameidx', simnameidx, ...
             'run', dump.run, ...
             'mean_cluster_speed', mean_cluster_speed, ...
             'mean_cluster_move', mean_cluster_move, ...
             'sd_cluster_move', sd_cluster_move, ...
             'sd_cluster_speed', sd_cluster_speed, ...
             'cluster_count', cluster_count, ...
             'avgcoverage', avgcoverage, ...
             'cumcoverage', cumcoverage, ...
             'mean_cluster_size', mean_cluster_size, ... 
             'sd_cluster_size', sd_cluster_size, ... 
             'mean_from_truth', mean_cluster_fromtruth, ... 
             'sd_from_truth', sd_cluster_fromtruth, ...
             'clusters_size', clusters_size, ...
             'clusters_fromtruth', clusters_fromtruth, ...
             'clusters_speed', clusters_speed, ...
             'clusters_move', clusters_move ...
        );
        
        
        
        if (PLOTS)
            
            close all;
            figure;
    
            subplot(2,3,1);
            hold on
            plot(1:nIter, avgcoverage, 'r')
            plot(1:nIter, cumcoverage, 'b')
            title('Avg and Cum coverage');
            hold off

            subplot(2,3,2);
            hold on
            plot(1:nIter, mean_cluster_speed, 'r')
            plot(1:nIter, sd_cluster_speed, 'b')
            title('Mean and Std cluster speed');
            hold off

            subplot(2,3,3);
            hold on
            plot(1:nIter, mean_cluster_move, 'r')
            plot(1:nIter, sd_cluster_move, 'b')
            title('Mean and Std cluster move');
            hold off

            subplot(2,3,4);
            hold on
            plot(1:nIter, cluster_count, 'g')
            plot(1:nIter, mean_cluster_size, 'r')
            plot(1:nIter, sd_cluster_size, 'b')
            title('Cluster Count, Mean and Std Size');
            hold off

            subplot(2,3,5);
            hold on
            plot(1:nIter, mean_cluster_fromtruth, 'r')
            plot(1:nIter, sd_cluster_fromtruth, 'b')
            title('Mean and Std cluster from truth');
            hold off

            
            % Creating a string with the description of the parameters

            paramString = format_sim_params_for_plot_display(simName, NAME, dump.parameters);

            annotation('textbox', [0.7, 0.45, 0, 0], 'string', paramString, ...
                'BackgroundColor', 'white', ...
                'EdgeColor', 'black', ...
                'LineStyle', '-' ...
            );

            waitforbuttonpress;
            
            figure;
            
            % Eliminating groups with only 1 agent
            % The first iteration (t0) is removed, there is no influence among
            % clusters 
            valid_clusters_idx = cell(nIter-1);
            
            subplot(1,3,1);
            hold on
            for z=2:nIter
                csizes = clusters_size{z};
                idxs = find(csizes > CLU_CUTOFF);
                valid_clusters_idx{z} = idxs;
                plot(clusters_size{z}(idxs), clusters_speed{z}(idxs), 'rx');
            end
            hold off
            title('Cluster size vs speed');
            
            subplot(1,3,2);
            hold on
            for z=2:nIter
                idxs = valid_clusters_idx{z};
                plot(clusters_size{z}(idxs), clusters_fromtruth{z}(idxs), 'rx');
            end
            hold off
            title('Cluster size vs fromtruth');

            subplot(1,3,3);
            hold on
            for z=2:nIter
                idxs = valid_clusters_idx{z};
                plot(clusters_size{z}(idxs), clusters_move{z}(idxs), 'rx');
            end
            hold off
            title('Cluster size vs move');
             
            waitforbuttonpress;

        end
        
        
        if (CSV_DUMP)
            
            simData = summaryObj{validFileIdx};
            
             % SAVING ONLY EVERY X ITERATIONS        
             for k = 1:size(idxsIters,2)
                z = idxsIters(k);
                stepData = simData(z);
                % 1 Line
                clu_macro_string = csv_format_row_clusters_macro(stepData, simName, z);
                fprintf(fidClustersMacro,'%s\n', clu_macro_string);
   
                % SAVING ONLY CLUSTERS of SIZE > CUTOFF
                idxs = find(stepData.clusters_size > CLU_CUTOFF);
                
                % Multiple Lines
                for i=1:length(idxs)
                    clu_micro_string = csv_format_row_clusters_micro(stepData, simName, z, idxs(i)); 
                    fprintf(fidClustersMicro,'%s\n', clu_micro_string);   
                end
 
             end
             
        end
        
    end % File loop
    
    save([ dumpDir 'sums' ], 'global_count_sum', 'global_count_sumsquared');
    
    
    N = validFileIdx; % last value
    df = N - 1;
    
    % Computing global stats  
    t_count_avg = global_count_sum / nFiles; 
    t_count_sd = sqrt(((global_count_sumsquared - ((global_count_sum).^2 / N))) / df);
    t_count_se = t_count_sd / sqrt(nFiles);  
    t_count_ci = t_count_se * tquant(CI_INT, nFiles-1);
    
    t_cover_avg = global_coverage_sum / nFiles; 
    t_cover_sd = sqrt(((global_coverage_sumsquared - ((global_coverage_sum).^2 / N))) / df);
    t_cover_se = t_cover_sd / sqrt(nFiles);  
    t_cover_ci = t_cover_se * tquant(CI_INT, nFiles-1);
    
    t_cover_cum_avg = global_coverage_cum_sum / nFiles; 
    t_cover_cum_sd = sqrt(((global_coverage_cum_sumsquared - ((global_coverage_cum_sum).^2 / N))) / df);
    t_cover_cum_se = t_cover_cum_sd / sqrt(nFiles);  
    t_cover_cum_ci = t_cover_cum_se * tquant(CI_INT, nFiles-1);
    
    t_speed_avg = global_speed_sum / nFiles; 
    t_speed_sd = sqrt(((global_speed_sumsquared - ((global_speed_sum).^2 / N))) / df);
    t_speed_se = t_speed_sd / sqrt(nFiles);  
    t_speed_ci = t_speed_se * tquant(CI_INT, nFiles-1);
    
    t_move_avg = global_move_sum / nFiles; 
    t_move_sd = sqrt(((global_move_sumsquared - ((global_move_sum).^2 / N))) / df);
    t_move_se = t_move_sd / sqrt(nFiles);  
    t_move_ci = t_move_se * tquant(CI_INT, nFiles-1);
    
    t_size_avg = global_size_sum / nFiles; 
    t_size_sd = sqrt(((global_size_sumsquared - ((global_size_sum).^2 / N))) / df);
    t_size_se = t_size_sd / sqrt(nFiles);  
    t_size_ci = t_size_se * tquant(CI_INT, nFiles-1);

    t_fromtruth_avg = global_fromtruth_sum / nFiles; 
    t_fromtruth_sd = sqrt(((global_fromtruth_sumsquared - ((global_fromtruth_sum).^2 / N))) / df);
    t_fromtruth_se = t_fromtruth_sd / sqrt(nFiles);  
    t_fromtruth_ci = t_fromtruth_se * tquant(CI_INT, nFiles-1);
 
    
    % Save summaryObj
    % save([dumpDir 'temporalysis.mat'], 'summaryObj');
    
    if (CSV_DUMP)
        
        fclose(fidParam);
        fclose(fidClustersMacro);
        fclose(fidClustersMicro);
        
        % SAVING ONLY EVERY X ITERATIONS        
        for k = 1:size(idxsIters,2)
            z = idxsIters(k);
            clu_macro_avg_string = sprintf('"%s",%u,%u,%u,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f', ...
                simName, simnameidx, dump.run, t_count_avg(z), t_count_sd(z), t_count_se(z), t_count_ci(z), ...
                t_cover_avg(z), t_cover_sd(z), t_cover_se(z), t_cover_ci(z), ...
                t_cover_cum_avg(z), t_cover_cum_sd(z), t_cover_cum_se(z), t_cover_cum_ci(z), ...
                t_speed_avg(z), t_speed_sd(z), t_speed_se(z), t_speed_ci(z), ...
                t_move_avg(z), t_move_sd(z), t_move_se(z), t_move_ci(z), ...
                t_size_avg(z),t_size_sd(z), t_size_se(z), t_size_ci(z), ...
                t_fromtruth_avg(z), t_fromtruth_sd(z), t_fromtruth_se(z), t_fromtruth_ci(z));
            
            fprintf(fidClustersMacroAvg,'%s\n', clu_macro_avg_string);   
         end

        fclose(fidClustersMacroAvg);
    end
    
    
    
end

%%




