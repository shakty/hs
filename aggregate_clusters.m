function aggregate_clusters(dumpDir, subDir, outDir)
    tic;
    
    % This function aggregates the params as well.
    
    % Creating outDir if not existing.
    if (exist(outDir, 'dir') == 0)
        mkdir(outDir);
    end
    
    % AGENTS
    
    dirPath = [dumpDir subDir 'clusters/'];
    
    % PreLoad one file (to create the variable nIter).
    fileName = 'sums_1-1.mat';    
    load([dirPath fileName]);
    nIter = length(global_speed_sum);
    
    files = dir(dirPath);
    fileIndex = find(~[files.isdir]);

    if (isempty(fileIndex))
        error('Invalid Directory Selected');
    end

    % Number of files
    nFiles = length(fileIndex);    
    
    % Reset file index for statistics;
    validFileIdx = 0;
    
    % Init g_global stats arrays.
    g_global_count_sum = zeros(nIter,1);
    g_global_count_sumsquared = zeros(nIter,1);
    g_global_maxsize_sum = zeros(nIter,1);
    g_global_maxsize_sumsquared = zeros(nIter,1);
    g_global_meansize_sum = zeros(nIter,1);
    g_global_meansize_sumsquared = zeros(nIter,1); 
    g_global_speed_sum = zeros(nIter,1);
    g_global_speed_sumsquared = zeros(nIter,1);
    g_global_movs_sum = zeros(nIter,1);
    g_global_movs_sumsquared = zeros(nIter,1);
    g_global_fromtruth_sum = zeros(nIter,1);
    g_global_fromtruth_sumsquared = zeros(nIter,1);
    g_global_bigcpdist_sum = zeros(nIter,1);
    g_global_bigcpdist_sumsquared = zeros(nIter,1);
    
    headers_clusters_macro = {
        'simname', ...
        'simcount', ...
        'run', ...
        't', ...
        'count', ...
        'size.max', ...
        'size.avg', ...
        'size.sd', ...
        'speed.avg', ...
        'speed.sd', ...
        'move.avg', ...
        'move.sd', ...
        'fromtruth.avg', ...
        'fromtruth.sd', ...
        'bigc.pdist.mean', ...
        'bigc.pdist.sd' ...   
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

    clustersMacroFileName = [outDir 'clusters_macro.csv'];
    % This function overwrites exiting files.
    write_csv_headers(clustersMacroFileName, headers_clusters_macro);   
    
    % Scan all files (also those we are interested in.
    for f = 1:nFiles

        append = files(fileIndex(f)).name;
        fileName = [dirPath append];

        % We load only sums_**.mat
        [PATH, NAME, EXT] = fileparts(append);
        if (~strcmpi(EXT,'.mat') ...
            || strfind(NAME, 'sums_') ~= 1)            
            continue;
        end

        % Extracting the part 1-1.mat from sums_1-1.mat
        % 6 = length('sums_') + 1;
        simnameidx = [NAME(6:length(NAME)) EXT];
        
        % Merge the corresponding CSV file.        
        clustersFileCSV = [dirPath 'clusters_macro_' simnameidx '.csv'];
        mergeCommand = sprintf('cat %s >> %s', ...
            clustersFileCSV, clustersMacroFileName);        
        system(mergeCommand);
        
        % Load file to compute round statistics.
        load(fileName);
        
        % Increment index of valid files.
        validFileIdx = validFileIdx + 1;
        
        % Updates stats arrays.
        g_global_count_sum = g_global_count_sum + global_count_sum;
        g_global_count_sumsquared = g_global_count_sumsquared + global_count_sumsquared;
        g_global_meansize_sum = g_global_meansize_sum + global_meansize_sum;
        g_global_meansize_sumsquared = g_global_meansize_sumsquared + global_meansize_sumsquared;
        g_global_maxsize_sum = g_global_maxsize_sum + global_maxsize_sum;
        g_global_maxsize_sumsquared = g_global_maxsize_sumsquared + global_maxsize_sumsquared;        
        g_global_speed_sum = g_global_speed_sum + global_speed_sum;
        g_global_speed_sumsquared = g_global_speed_sumsquared + global_speed_sumsquared;
        g_global_movs_sum = g_global_movs_sum + global_movs_sum;
        g_global_movs_sumsquared = g_global_movs_sumsquared + global_movs_sumsquared;
        g_global_fromtruth_sum = g_global_fromtruth_sum + global_fromtruth_sum;
        g_global_fromtruth_sumsquared = g_global_fromtruth_sumsquared + global_fromtruth_sumsquared;
        g_global_bigcpdist_sum = g_global_bigcpdist_sum + global_bigc_pairwise_pdist_sum;
        g_global_bigcpdist_sumsquared = g_global_bigcpdist_sumsquared + global_bigc_pairwise_pdist_sumsquared;
        
        % Delete file variables to be sure we do not overwrite some by
        % mistake.
        clearvars global_*;
        
    end
    
    % N: observations.
    N = validFileIdx;
    % df: degree of freedom.
    df = N - 1;
    % 95 percent confidence interval.
    CI_INT = 0.95/2 + 0.5;
    
    % Computing g_global stats.
    
    t_count_avg = g_global_count_sum / N;
    t_count_sd = sqrt(((g_global_count_sumsquared - ((g_global_count_sum).^2 / N))) / df);
    t_count_se = t_count_sd / sqrt(N);
    t_count_ci = t_count_se * tquant(CI_INT, df);
    
    t_meansize_avg = g_global_meansize_sum / N; 
    t_meansize_sd = sqrt(((g_global_meansize_sumsquared - ((g_global_meansize_sum).^2 / N))) / df);
    t_meansize_se = t_meansize_sd / sqrt(N);  
    t_meansize_ci = t_meansize_se * tquant(CI_INT, df);
    
    t_maxsize_avg = g_global_maxsize_sum / N; 
    t_maxsize_sd = sqrt(((g_global_maxsize_sumsquared - ((g_global_maxsize_sum).^2 / N))) / df);
    t_maxsize_se = t_meansize_sd / sqrt(N);  
    t_maxsize_ci = t_meansize_se * tquant(CI_INT, df);
    
    t_speed_avg = g_global_speed_sum / N; 
    t_speed_sd = sqrt(((g_global_speed_sumsquared - ((g_global_speed_sum).^2 / N))) / df);
    t_speed_se = t_speed_sd / sqrt(N);  
    t_speed_ci = t_speed_se * tquant(CI_INT, df);
    
    t_move_avg = g_global_movs_sum / N; 
    t_move_sd = sqrt(((g_global_movs_sumsquared - ((g_global_movs_sum).^2 / N))) / df);
    t_move_se = t_move_sd / sqrt(N);  
    t_move_ci = t_move_se * tquant(CI_INT, df);

    t_fromtruth_avg = g_global_fromtruth_sum / N; 
    t_fromtruth_sd = sqrt(((g_global_fromtruth_sumsquared - ((g_global_fromtruth_sum).^2 / N))) / df);
    t_fromtruth_se = t_fromtruth_sd / sqrt(N);  
    t_fromtruth_ci = t_fromtruth_se * tquant(CI_INT, df);
    
    t_bigcpdist_avg = g_global_bigcpdist_sum / N; 
    t_bigcpdist_sd = sqrt(((g_global_bigcpdist_sumsquared - ((g_global_bigcpdist_sum).^2 / N))) / df);
    t_bigcpdist_se = t_move_sd / sqrt(N);
    t_bigcpdist_ci = t_move_se * tquant(CI_INT, df);
    
    headers_clusters_macro_avg = {
        'simname', ...
        'simcount', ...
        't', ...
        'count.avg', ...
        'count.sd', ...
        'count.se', ...
        'count.ci', ...
        'meansize.avg', ...
        'meansize.sd', ...
        'meansize.se', ...
        'meansize.ci', ...
        'maxsize.avg', ...
        'maxsize.sd', ...
        'maxsize.se', ...
        'maxsize.ci', ...
        'speed.avg', ...
        'speed.sd', ...
        'speed.se', ...
        'speed.ci', ...
        'move.avg', ...
        'move.sd', ...
        'move.se', ...
        'move.ci', ...
        'fromtruth.avg', ...
        'fromtruth.sd', ...
        'fromtruth.se', ...
        'fromtruth.ci', ...
        'bigc.pdist.avg', ...
        'bigc.pdist.sd', ...
        'bigc.pdist.se', ...
        'bigc.pdist.ci' ...
    };
    
    clustersAvgFileName = [outDir 'clusters_avg_all.csv'];
    write_csv_headers(clustersAvgFileName, headers_clusters_macro_avg);
    fidClustersMacroAvg = fopen(clustersAvgFileName, 'a');
     
    for z = 1:nIter
        clu_macro_avg_string = sprintf('"%s",%u,%u,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f', ...
            subDir, N, z, ...
            t_count_avg(z), t_count_sd(z), t_count_se(z), t_count_ci(z), ...
            t_meansize_avg(z),t_meansize_sd(z), t_meansize_se(z), t_meansize_ci(z), ...
            t_maxsize_avg(z),t_maxsize_sd(z), t_maxsize_se(z), t_maxsize_ci(z), ...
            t_speed_avg(z), t_speed_sd(z), t_speed_se(z), t_speed_ci(z), ...
            t_move_avg(z), t_move_sd(z), t_move_se(z), t_move_ci(z), ...
            t_fromtruth_avg(z), t_fromtruth_sd(z), t_fromtruth_se(z), t_fromtruth_ci(z), ...
            t_bigcpdist_avg(z),t_bigcpdist_sd(z), t_bigcpdist_se(z), t_bigcpdist_ci(z));
        
        fprintf(fidClustersMacroAvg, '%s\n', clu_macro_avg_string);   
     end

    fclose(fidClustersMacroAvg);
    
    % Delete the files ?
end
