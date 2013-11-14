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
        'attrtype', ...
        'attr_on_v', ...
        'seed', ...
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
    % idxsIters = find(mod(1:nIter, DUMP_RATE) == 0);
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
        [PATH,NAME_tmp,EXT_tmp] = fileparts(fileName);
        if (~strcmpi(EXT_tmp,'.mat') ...
            || strcmp(NAME_tmp, 'sums') == 1 ...   
            || strcmp(NAME_tmp, 'temporalysis') == 1 ...
            || strcmp(NAME_tmp, 'global_count_sum') == 1 ...
            || strcmp(NAME_tmp, 'global_count_sumsquared') == 1 )
            continue;
        end

        % It was a valid file, so update the NAME and EXT var
        NAME = NAME_tmp;
        EXT = EXT_tmp;
        
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
        
        % vector of movements of agents at time t0
        movs = zeros(size(pos,2), 1);
        mean_cluster_move = 0;
        sd_cluster_move = 0;
        
        for i = 1:nIter
            %avg speed
            mean_cluster_speed = mean(colnorm(v(:,:,i),2));
            % avg movement
            if (i>1)
                movs = abs(pos(:,:,i)-pos(:,:,i-1));
                mean_cluster_move = mean(mean(movs));
            end
            
            %avg share of space
            coverage_matrix = countAgents(pos(:,:,i), PRECISION);
            avgcoverage = nnz(coverage_matrix) / TOTAL_CELLS;

            %cum share of space
            if (i==1)
                cum_coverage_matrix = coverage_matrix;
            else
                cum_coverage_matrix= coverage_matrix + cum_coverage_matrix;
            end
            cumcoverage = nnz(cum_coverage_matrix) / TOTAL_CELLS;

            try
                % Z the results of HCLUST 
                % T the cluster of each agent
                % C the number of clusters            
                [Z, T, C] = clusterize(pos(:,:,i));
                cluster_count = C;

                [d, Gc, AvgGDist, avgGroupSpeed, avgGroupMove] = cluster_stats(T, dump.truth, pos(:,:,i), v(:,:,i), movs);

                mean_cluster_size = mean(Gc);
                sd_cluster_size = std(Gc);

                mean_cluster_fromtruth = mean(AvgGDist);
                sd_cluster_fromtruth = std(AvgGDist);

                sd_cluster_speed = std(avgGroupSpeed);
                sd_cluster_move = std(avgGroupMove);

                clusters_size = Gc;
                clusters_fromtruth = AvgGDist;
                clusters_speed = avgGroupSpeed;
                clusters_move = avgGroupMove;

                % SUMMING UP AVG statistics
                global_count_sum(i) = global_count_sum(i) + C;
                global_count_sumsquared(i) = global_count_sumsquared(i) + C^2;
                global_coverage_sum(i) = global_coverage_sum(i) + avgcoverage;
                global_coverage_sumsquared(i) = global_coverage_sumsquared(i) + avgcoverage^2;
                global_coverage_cum_sum(i) = global_coverage_cum_sum(i) + cumcoverage;
                global_coverage_cum_sumsquared(i) = global_coverage_cum_sumsquared(i) + cumcoverage^2;
                global_speed_sum(i) = global_speed_sum(i) + mean_cluster_speed;
                global_speed_sumsquared(i) = global_speed_sumsquared(i) + mean_cluster_speed^2;
                global_move_sum(i) = global_move_sum(i) + mean_cluster_move;
                global_move_sumsquared(i) = global_move_sumsquared(i) + mean_cluster_move^2;
                global_size_sum(i) = global_size_sum(i) + mean_cluster_size;
                global_size_sumsquared(i) = global_size_sumsquared(i) + mean_cluster_size^2;
                global_fromtruth_sum(i) = global_fromtruth_sum(i) +  mean_cluster_fromtruth;
                global_fromtruth_sumsquared(i) = global_fromtruth_sumsquared(i) + mean_cluster_fromtruth^2;
                
            catch err
                
                err
                sprintf('Error at iter: %i, fileIdx: %i, fileName: %s', i, validFileIdx, fileName)
                
                mean_cluster_size = 0;
                sd_cluster_size = 0;

                mean_cluster_fromtruth = 0;
                sd_cluster_fromtruth = 0;

                sd_cluster_speed = 0;
                sd_cluster_move = 0;

                clusters_size = 0;
                clusters_fromtruth = 0;
                clusters_speed = 0;
                clusters_move = 0;
            
            end
            
         if (CSV_DUMP)
            
            
             % SAVING ONLY EVERY X ITERATIONS        
             if (mod(i, DUMP_RATE) == 0)
                                
                 stepData = struct(...
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
            
            

        end
        
% Trying to save memory. We dump in real time.        
%         summaryObj{validFileIdx} = struct(...
%              'simnameidx', simnameidx, ...
%              'run', dump.run, ...
%              'mean_cluster_speed', mean_cluster_speed, ...
%              'mean_cluster_move', mean_cluster_move, ...
%              'sd_cluster_move', sd_cluster_move, ...
%              'sd_cluster_speed', sd_cluster_speed, ...
%              'cluster_count', cluster_count, ...
%              'avgcoverage', avgcoverage, ...
%              'cumcoverage', cumcoverage, ...
%              'mean_cluster_size', mean_cluster_size, ... 
%              'sd_cluster_size', sd_cluster_size, ... 
%              'mean_from_truth', mean_cluster_fromtruth, ... 
%              'sd_from_truth', sd_cluster_fromtruth, ...
%              'clusters_size', clusters_size, ...
%              'clusters_fromtruth', clusters_fromtruth, ...
%              'clusters_speed', clusters_speed, ...
%              'clusters_move', clusters_move ...
%         );   
        
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
        
% Trying to use less memory. We dump in real time, and not at the end       
%         if (CSV_DUMP)
%             
%             simData = summaryObj{validFileIdx};
%             
%              % SAVING ONLY EVERY X ITERATIONS        
%              for k = 1:size(idxsIters,2)
%                 z = idxsIters(k);
%                 stepData = simData(z);
%                 % 1 Line
%                 clu_macro_string = csv_format_row_clusters_macro(stepData, simName, z);
%                 fprintf(fidClustersMacro,'%s\n', clu_macro_string);
%    
%                 % SAVING ONLY CLUSTERS of SIZE > CUTOFF
%                 idxs = find(stepData.clusters_size > CLU_CUTOFF);
%                 
%                 % Multiple Lines
%                 for i=1:length(idxs)
%                     clu_micro_string = csv_format_row_clusters_micro(stepData, simName, z, idxs(i)); 
%                     fprintf(fidClustersMicro,'%s\n', clu_micro_string);   
%                 end
%  
%              end
%              
%         end
        
    end % File loop
    
    N = validFileIdx; % last value
    df = N - 1;
    
    save([ dumpDir 'sums' ], 'global_count_sum', ...
                             'global_count_sumsquared', ...
                             'global_coverage_sum', ...
                             'global_coverage_sumsquared', ...
                             'global_coverage_cum_sum', ...
                             'global_coverage_cum_sumsquared', ...
                             'global_speed_sum', ...
                             'global_speed_sumsquared', ...
                             'global_move_sum', ...
                             'global_move_sumsquared', ...
                             'global_size_sum', ...
                             'global_size_sumsquared', ...
                             'global_fromtruth_sum', ...
                             'global_fromtruth_sumsquared', ...
                             'N' ... % number of files
    );
    
    
    % Computing global stats  
    t_count_avg = global_count_sum / N; 
    t_count_sd = sqrt(((global_count_sumsquared - ((global_count_sum).^2 / N))) / df);
    t_count_se = t_count_sd / sqrt(N);  
    t_count_ci = t_count_se * tquant(CI_INT, df);
    
    t_cover_avg = global_coverage_sum / N; 
    t_cover_sd = sqrt(((global_coverage_sumsquared - ((global_coverage_sum).^2 / N))) / df);
    t_cover_se = t_cover_sd / sqrt(N);  
    t_cover_ci = t_cover_se * tquant(CI_INT, df);
    
    t_cover_cum_avg = global_coverage_cum_sum / N; 
    t_cover_cum_sd = sqrt(((global_coverage_cum_sumsquared - ((global_coverage_cum_sum).^2 / N))) / df);
    t_cover_cum_se = t_cover_cum_sd / sqrt(N);  
    t_cover_cum_ci = t_cover_cum_se * tquant(CI_INT, df);
    
    t_speed_avg = global_speed_sum / N; 
    t_speed_sd = sqrt(((global_speed_sumsquared - ((global_speed_sum).^2 / N))) / df);
    t_speed_se = t_speed_sd / sqrt(N);  
    t_speed_ci = t_speed_se * tquant(CI_INT, df);
    
    t_move_avg = global_move_sum / N; 
    t_move_sd = sqrt(((global_move_sumsquared - ((global_move_sum).^2 / N))) / df);
    t_move_se = t_move_sd / sqrt(N);  
    t_move_ci = t_move_se * tquant(CI_INT, df);
    
    t_size_avg = global_size_sum / N; 
    t_size_sd = sqrt(((global_size_sumsquared - ((global_size_sum).^2 / N))) / df);
    t_size_se = t_size_sd / sqrt(N);  
    t_size_ci = t_size_se * tquant(CI_INT, df);

    t_fromtruth_avg = global_fromtruth_sum / N; 
    t_fromtruth_sd = sqrt(((global_fromtruth_sumsquared - ((global_fromtruth_sum).^2 / N))) / df);
    t_fromtruth_se = t_fromtruth_sd / sqrt(N);  
    t_fromtruth_ci = t_fromtruth_se * tquant(CI_INT, df);
 
    
    % Save summaryObj
    % save([dumpDir 'temporalysis.mat'], 'summaryObj');
    
    if (CSV_DUMP)
        
        fclose(fidParam);
        fclose(fidClustersMacro);
        fclose(fidClustersMicro);
        
        % SAVING ALL ITERATIONS for the AVG
        for z = 1:nIter
            clu_macro_avg_string = sprintf('"%s",%u,%u,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f', ...
                simName, N, z, t_count_avg(z), t_count_sd(z), t_count_se(z), t_count_ci(z), ...
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




