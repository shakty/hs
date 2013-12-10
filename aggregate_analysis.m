function aggregate_analysis(dumpDir, subDir, outDir)
    tic;
    
%     headers_clusters_macro_avg = {
%         'simname', ...
%         'simcount', ...
%         't', ...
%         'count.avg', ...
%         'count.sd', ...
%         'count.se', ...
%         'count.ci', ...
%         'coverage.avg', ...
%         'coverage.sd', ...
%         'coverage.se', ...
%         'coverage.ci', ...
%         'coverage.cum.avg', ...
%         'coverage.cum.sd', ...
%         'coverage.cum.se', ...
%         'coverage.cum.ci', ...
%         'speed.avg', ...
%         'speed.sd', ...
%         'speed.se', ...
%         'speed.ci', ...
%         'move.avg', ...
%         'move.sd', ...
%         'move.se', ...
%         'move.ci', ...
%         'size.avg', ...
%         'size.sd', ...
%         'size.se', ...
%         'size.ci', ...
%         'fromtruth.avg', ...
%         'fromtruth.sd', ...
%         'fromtruth.se', ...
%         'fromtruth.ci' ...
%     };
% 
%     CI_INT = 0.95/2 + 0.5;
% 
%     % Creating variables used to save results to file.
%     
%     headers_truthradius_avg = {
%         'simname', ...
%         'simcount', ...
%         't' ...
%     };
% 
%     avg_symbols = { '_mean', '_sd', '_se', '_ci' };
%     nSyms = length(avg_symbols);
% 
%     % nRadiuses
%     nRadiuses = length(RADIUSs);
%     hAvgStartFrom = 3;
%     radiusesStr = {};
%     for i = 1 : nRadiuses
%         radiusStr = ['r_' num2str(RADIUSs(i))];
%         radiusesStr{i} = radiusStr;
%         for j = 1 : nSyms
%             idxJ = (i -1) * nSyms + j + hAvgStartFrom;
%             radiusesStrAvg = [radiusStr avg_symbols{j}];
%             headers_truthradius_avg{idxJ} = radiusesStrAvg;
%         end
%         
%     end
%     
%     % Adding not in radius
%     i = i + 1;
%     radiusStr = 'r_out';
%     radiusesStr{i} = radiusStr;
%     
%     for j = 1 : nSyms
%         idxJ = (i -1) * nSyms + j + hAvgStartFrom;
%         headers_truthradius_avg{idxJ} = [radiusStr avg_symbols{j}];
%     end
%     
%     % Adding flag consensus on Truth.
%     i = i + 1;
%     radiusStr = 'consensus';
%     radiusesStr{i} = radiusStr;
%     
%     for j = 1 : nSyms
%         idxJ = (i -1) * nSyms + j + hAvgStartFrom;
%         headers_truthradius_avg{idxJ} = [radiusStr avg_symbols{j}];
%     end
%     
%     % +1 is NOT_IN_RADIUS.
%     nRadiusesPlusOne = nRadiuses + 1;
%     
%     % End creating variables used to save results to file.
%     
%     
%     % GLOBAL variables
%   
%     globalRadiusCounts = zeros(nIter, nRadiusesPlusOne);
%     globalRadiusCounts_squared = zeros(nIter, nRadiusesPlusOne);
%     globalConsensusOnTruth = zeros(1, nIter);
%     globalConsensusOnTruth_squared = zeros(1, nIter);
%     
%     wait(j);
%     
%     % Delete Job
%     delete(j);
%     
%     
%     files = dir(outDir);
%     fileIndex = find(~[files.isdir]);
% 
%     if (isempty(fileIndex))
%         error('Invalid Directory Selected');
%     end
%     
%     nFiles = length(fileIndex);
    
    % PARAMS (for all)
    
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
    
    paramsFileName = [outDir 'params.csv'];
    % This function overwrites exiting files.
    write_csv_headers(paramsFileName, headers_params);
    
    % AGENTS
    
    % Reset file index for statistics;
    validFileIdx = 0;
    
    % Init global stats arrays.
    g_global_coverage_sum = zeros(nIter,1);
    g_global_coverage_sumsquared = zeros(nIter,1);
    g_global_coverage_cum_sum = zeros(nIter,1);
    g_global_coverage_cum_sumsquared = zeros(nIter,1);
    g_global_speed_sum = zeros(nIter,1);
    g_global_speed_sumsquared = zeros(nIter,1);
    g_global_movs_sum = zeros(nIter,1);
    g_global_movs_sumsquared = zeros(nIter,1); 
    g_global_fromtruth_sum = zeros(nIter,1);
    g_global_fromtruth_sumsquared = zeros(nIter,1);
    g_global_pdist_sum =  zeros(nIter,1);
    g_global_pdist_sumsquared =  zeros(nIter,1);
    
    
    headers_agents = {
        'simname', ...
        'simcount', ...
        'run', ...
        't', ...
        'coverage', ...
        'coverage.cum', ...
        'speed.avg', ...
        'speed.sd', ...
        'move.avg', ...
        'move.sd', ...
        'fromtruth.avg', ...
        'fromtruth.sd', ...
        'pdist.mean', ...
        'pdist.sd' ...   
    };

    agentsFileName = [outDir 'agents.csv'];
    % This function overwrites exiting files.
    write_csv_headers(agentsFileName, headers_agents);
    
    
    % Scan all files (also those we are interested in.
    for f = 1:nFiles

        append = files(fileIndex(f)).name;
        fileName = [dumpDir, append];

        % We load only .mat
        [PATH, NAME, EXT] = fileparts(fileName);
        if (~strcmpi(EXT,'.mat') ...
            || strcmp(NAME, 'sums') == 1 ...   
            || strcmp(NAME, 'temporalysis') == 1 ...
            || strcmp(NAME, 'global_count_sum') == 1 ...
            || strcmp(NAME, 'global_count_sumsquared') == 1 )
            continue;
        end

        % Merge the corresponding CSV file.        
        agentsFileCSV = ['agents_' NAME EXT '.csv'];
        mergeCommand = sprintf('%s >> %s', agentsFileCSV, agentsFileName);        
        system(mergeCommand);
        
        % Just in this case, merge also the params file.
        paramsFileCSV = ['params_' NAME EXT '.csv'];
        mergeCommand = sprintf('%s >> %s', paramsFileCSV, paramsFileName);        
        system(mergeCommand);
        
        % Load file to compute round statistics.
        load(fileName);
        
        % Increment index of valid files.
        validFileIdx = validFileIdx + 1;
        
        % Updates stats arrays.
        g_global_count_sum = g_global_count_sum + global_count_sum;
        g_global_count_sumsquared = g_global_count_sumsquared + global_count_sumsquared;
        g_global_coverage_sum = g_global_coverage_sum + global_coverage_sum;
        g_global_coverage_sumsquared = g_global_coverage_sumsquared + global_coverage_sumsquared;
        g_global_coverage_cum_sum = g_global_coverage_cum_sum + global_coverage_cum_sum;
        g_global_coverage_cum_sumsquared = g_global_coverage_cum_sumsquared + global_coverage_cum_sumsquared;
        g_global_speed_sum = g_global_speed_sum + global_speed_sum;
        g_global_speed_sumsquared = g_global_speed_sumsquared + global_speed_sumsquared;
        g_global_movs_sum = g_global_movs_sum + global_move_sum;
        g_global_movs_sumsquared = g_global_movs_sumsquared + global_move_sumsquared;
        g_global_size_sum = g_global_size_sum + global_size_sum;
        g_global_size_sumsquared = g_global_size_sumsquared + global_size_sumsquared;
        g_global_fromtruth_sum = g_global_fromtruth_sum + global_fromtruth_sum;
        g_global_fromtruth_sumsquared = g_global_fromtruth_sumsquared + global_fromtruth_sumsquared;     
        g_global_pdist_sum = g_global_pdist_sum + global_size_sum;
        g_global_pdist_sumsquared = g_global_pdist_sumsquared + global_size_sumsquared;
    end
    
    % N: observations.
    N = validFileIdx;
    % df: degree of freedom.
    df = N - 1;
    
    % Computing g_global stats.
    t_count_avg = g_global_count_sum / N; 
    t_count_sd = sqrt(((g_global_count_sumsquared - ((g_global_count_sum).^2 / N))) / df);
    t_count_se = t_count_sd / sqrt(N);  
    t_count_ci = t_count_se * tquant(CI_INT, df);
    
    t_cover_avg = g_global_coverage_sum / N; 
    t_cover_sd = sqrt(((g_global_coverage_sumsquared - ((g_global_coverage_sum).^2 / N))) / df);
    t_cover_se = t_cover_sd / sqrt(N);  
    t_cover_ci = t_cover_se * tquant(CI_INT, df);
    
    t_cover_cum_avg = g_global_coverage_cum_sum / N; 
    t_cover_cum_sd = sqrt(((g_global_coverage_cum_sumsquared - ((g_global_coverage_cum_sum).^2 / N))) / df);
    t_cover_cum_se = t_cover_cum_sd / sqrt(N);  
    t_cover_cum_ci = t_cover_cum_se * tquant(CI_INT, df);
    
    t_speed_avg = g_global_speed_sum / N; 
    t_speed_sd = sqrt(((g_global_speed_sumsquared - ((g_global_speed_sum).^2 / N))) / df);
    t_speed_se = t_speed_sd / sqrt(N);  
    t_speed_ci = t_speed_se * tquant(CI_INT, df);
    
    t_move_avg = g_global_move_sum / N; 
    t_move_sd = sqrt(((g_global_move_sumsquared - ((g_global_move_sum).^2 / N))) / df);
    t_move_se = t_move_sd / sqrt(N);  
    t_move_ci = t_move_se * tquant(CI_INT, df);
    
    t_size_avg = g_global_size_sum / N; 
    t_size_sd = sqrt(((g_global_size_sumsquared - ((g_global_size_sum).^2 / N))) / df);
    t_size_se = t_size_sd / sqrt(N);  
    t_size_ci = t_size_se * tquant(CI_INT, df);

    t_fromtruth_avg = g_global_fromtruth_sum / N; 
    t_fromtruth_sd = sqrt(((g_global_fromtruth_sumsquared - ((g_global_fromtruth_sum).^2 / N))) / df);
    t_fromtruth_se = t_fromtruth_sd / sqrt(N);  
    t_fromtruth_ci = t_fromtruth_se * tquant(CI_INT, df);
    
    t_pdist_avg = g_global_pdist_sum / N; 
    t_pdist_sd = sqrt(((g_global_pdist_sumsquared - ((g_global_pdist_sum).^2 / N))) / df);
    t_pdist_se = t_pdist_sd / sqrt(N);  
    t_pdist_ci = t_pdist_se * tquant(CI_INT, df);
    
    agentsAvgFileName = [outDir 'agents_avg_all.csv'];
    write_csv_headers(agentsAvgFileName, headers_agents_avg);
    fidAgentsAvg = fopen(agentsAvgFileName, 'a');
    
    % Saving all iterations
    for z = 1:nIter
        agents_avg_string = sprintf('"%s",%u,%u,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f%.4f,%.4f,%.4f,%.4f', ...
            subDir, totalFiles, z, t_count_avg(z), t_count_sd(z), t_count_se(z), t_count_ci(z), ...
            t_cover_avg(z), t_cover_sd(z), t_cover_se(z), t_cover_ci(z), ...
            t_cover_cum_avg(z), t_cover_cum_sd(z), t_cover_cum_se(z), t_cover_cum_ci(z), ...
            t_speed_avg(z), t_speed_sd(z), t_speed_se(z), t_speed_ci(z), ...
            t_move_avg(z), t_move_sd(z), t_move_se(z), t_move_ci(z), ...
            t_size_avg(z),t_size_sd(z), t_size_se(z), t_size_ci(z), ...
            t_fromtruth_avg(z), t_fromtruth_sd(z), t_fromtruth_se(z), t_fromtruth_ci(z), ...
            t_pdist_avg(z), t_pdist_sd(z), t_pdist_se(z), t_pdist_ci(z) ...
        );

        fprintf(fidAgentsAvg,'%s\n', agents_avg_string);   
     end

    fclose(fidAgentsAvg);
    
    % Delete the files ?    
    
    % RADIUS
    
    % Computing global stats            
    for i = 1 : nRadiusesPlusOne
        
            v = genvarname(['t_' radiusesStr(i) '_mean']);
            v = globalRadiusCounts(i) / N; 
            t_cover_sd = sqrt(((global_coverage_sumsquared - ((global_radiusCounts).^2 / N))) / df);
            t_cover_se = t_cover_sd / sqrt(N);  
            t_cover_ci = t_cover_se * tquant(CI_INT, df);
        
    end
 
    if (CSV_DUMP)
        
        fclose(fidTruthRadius);
        
        % SAVING ALL ITERATIONS for the AVG
        for z = 1:nIter
            % 23 elements
            agents_avg_string = sprintf('"%s",%u,%u,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f', ...
                simName, N, z, ...
                t_cover_avg(z), t_cover_sd(z), t_cover_se(z), t_cover_ci(z), ...
                t_cover_cum_avg(z), t_cover_cum_sd(z), t_cover_cum_se(z), t_cover_cum_ci(z), ...
                t_speed_avg(z), t_speed_sd(z), t_speed_se(z), t_speed_ci(z), ...
                t_movs_avg(z), t_movs_sd(z), t_movs_se(z), t_movs_ci(z), ...
                t_fromtruth_avg(z), t_fromtruth_sd(z), t_fromtruth_se(z), t_fromtruth_ci(z));
            
            fprintf(fidTruthRadiusAvg,'%s\n', agents_avg_string);   
         end

        fclose(fidTruthRadiusAvg);
    end
    
    % Delete the directory outDir

end
  
% 
%    
%     
% 
%   
%         % Aggregating Radius
%         
%         load([ dirPath '/radius/radius_' fileName '.csv']);
%         
%         
%         
%         % Aggregating Agents
%         
%         
%         % Aggregating Clusters
%         
%         load([ dirPath '/' 'sums' ]);
%         
%         totalFiles = totalFiles + N; % is found in sums
%         
%      
%         
%         g_global_count_sum = g_global_count_sum + global_count_sum;
%         g_global_count_sumsquared = g_global_count_sumsquared + global_count_sumsquared;
%         g_global_coverage_sum = g_global_coverage_sum + global_coverage_sum;
%         g_global_coverage_sumsquared = g_global_coverage_sumsquared + global_coverage_sumsquared;
%         g_global_coverage_cum_sum = g_global_coverage_cum_sum + global_coverage_cum_sum;
%         g_global_coverage_cum_sumsquared = g_global_coverage_cum_sumsquared + global_coverage_cum_sumsquared;
%         g_global_speed_sum = g_global_speed_sum + global_speed_sum;
%         g_global_speed_sumsquared = g_global_speed_sumsquared + global_speed_sumsquared;
%         g_global_move_sum = g_global_move_sum + global_move_sum;
%         g_global_move_sumsquared = g_global_move_sumsquared + global_move_sumsquared;
%         g_global_size_sum = g_global_size_sum + global_size_sum;
%         g_global_size_sumsquared = g_global_size_sumsquared + global_size_sumsquared;
%         g_global_fromtruth_sum = g_global_fromtruth_sum + global_fromtruth_sum;
%         g_global_fromtruth_sumsquared = g_global_fromtruth_sumsquared + global_fromtruth_sumsquared;     
% 
%     end
% 
% 
%     
%     
%     
%   %% Clusters
%   
%      
%     % Computing global stats  
%     t_count_avg = global_count_sum / N; 
%     t_count_sd = sqrt(((global_count_sumsquared - ((global_count_sum).^2 / N))) / df);
%     t_count_se = t_count_sd / sqrt(N);  
%     t_count_ci = t_count_se * tquant(CI_INT, df);
%     
%     t_cover_avg = global_coverage_sum / N; 
%     t_cover_sd = sqrt(((global_coverage_sumsquared - ((global_coverage_sum).^2 / N))) / df);
%     t_cover_se = t_cover_sd / sqrt(N);  
%     t_cover_ci = t_cover_se * tquant(CI_INT, df);
%     
%     t_cover_cum_avg = global_coverage_cum_sum / N; 
%     t_cover_cum_sd = sqrt(((global_coverage_cum_sumsquared - ((global_coverage_cum_sum).^2 / N))) / df);
%     t_cover_cum_se = t_cover_cum_sd / sqrt(N);  
%     t_cover_cum_ci = t_cover_cum_se * tquant(CI_INT, df);
%     
%     t_speed_avg = global_speed_sum / N; 
%     t_speed_sd = sqrt(((global_speed_sumsquared - ((global_speed_sum).^2 / N))) / df);
%     t_speed_se = t_speed_sd / sqrt(N);  
%     t_speed_ci = t_speed_se * tquant(CI_INT, df);
%     
%     t_move_avg = global_move_sum / N; 
%     t_move_sd = sqrt(((global_move_sumsquared - ((global_move_sum).^2 / N))) / df);
%     t_move_se = t_move_sd / sqrt(N);  
%     t_move_ci = t_move_se * tquant(CI_INT, df);
%     
%     t_size_avg = global_size_sum / N; 
%     t_size_sd = sqrt(((global_size_sumsquared - ((global_size_sum).^2 / N))) / df);
%     t_size_se = t_size_sd / sqrt(N);  
%     t_size_ci = t_size_se * tquant(CI_INT, df);
% 
%     t_fromtruth_avg = global_fromtruth_sum / N; 
%     t_fromtruth_sd = sqrt(((global_fromtruth_sumsquared - ((global_fromtruth_sum).^2 / N))) / df);
%     t_fromtruth_se = t_fromtruth_sd / sqrt(N);  
%     t_fromtruth_ci = t_fromtruth_se * tquant(CI_INT, df);
%  
%     
%     % Save summaryObj
%     % save([dumpDir 'temporalysis.mat'], 'summaryObj');
%     
%     if (CSV_DUMP)
%         
%        
%         fclose(fidClustersMacro);
%         fclose(fidClustersMicro);
%         
%         % SAVING ALL ITERATIONS for the AVG
%         for z = 1:nIter
%             clu_macro_avg_string = sprintf('"%s",%u,%u,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f', ...
%                 simName, N, z, t_count_avg(z), t_count_sd(z), t_count_se(z), t_count_ci(z), ...
%                 t_cover_avg(z), t_cover_sd(z), t_cover_se(z), t_cover_ci(z), ...
%                 t_cover_cum_avg(z), t_cover_cum_sd(z), t_cover_cum_se(z), t_cover_cum_ci(z), ...
%                 t_speed_avg(z), t_speed_sd(z), t_speed_se(z), t_speed_ci(z), ...
%                 t_move_avg(z), t_move_sd(z), t_move_se(z), t_move_ci(z), ...
%                 t_size_avg(z),t_size_sd(z), t_size_se(z), t_size_ci(z), ...
%                 t_fromtruth_avg(z), t_fromtruth_sd(z), t_fromtruth_se(z), t_fromtruth_ci(z));
%             
%             fprintf(fidClustersMacroAvg,'%s\n', clu_macro_avg_string);   
%          end
% 
%         fclose(fidClustersMacroAvg);
%     end
  
    