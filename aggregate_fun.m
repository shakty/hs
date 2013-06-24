function [] = temporalysis_fun(dumpDir, subDir)

    %% Saves simulations into properly formatted CSV files

    %% Add other directories to path
    path(path,'util/'); % Help functions
    path(path,'lib/'); % Help functions

    % Change default axes fonts.
    set(0,'DefaultAxesFontName', 'Times New Roman')
    set(0,'DefaultAxesFontSize', 14)

    % Change default text fonts.
    set(0,'DefaultTextFontname', 'Times New Roman')
    set(0,'DefaultTextFontSize', 14)

    
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
    

    DUMPDIR = [dumpDir subDir];
    
    dirs = dir(DUMPDIR);
    dirIndex = find([dirs.isdir]);

    if (isempty(dirIndex))
        error('Invalid Directory Selected');
    end
    
    totalFiles = 0;
    
    validFiles = 0;
    % Load all parameters matrices in one
    for d = 1:length(dirIndex)
        
        myDir = dirs(dirIndex(d)).name;
        
        if (strcmpi(myDir,'.') || strcmpi(myDir,'..'))
            continue;
        end
        
        validFiles = validFiles + 1;
        
        dirPath = [DUMPDIR myDir];

        load([ dirPath '/' 'sums' ]);
        
        totalFiles = totalFiles + N; % is found in sums
        
        % If it is the first open file, initialize the arrays
        if (validFiles == 1)
            nIter = length(global_count_sum);
            g_global_count_sum = zeros(nIter,1);
            g_global_count_sumsquared = zeros(nIter,1);
            g_global_coverage_sum = zeros(nIter,1);
            g_global_coverage_sumsquared = zeros(nIter,1);
            g_global_coverage_cum_sum = zeros(nIter,1);
            g_global_coverage_cum_sumsquared = zeros(nIter,1);
            g_global_speed_sum = zeros(nIter,1);
            g_global_speed_sumsquared = zeros(nIter,1);
            g_global_move_sum = zeros(nIter,1);
            g_global_move_sumsquared = zeros(nIter,1);
            g_global_size_sum = zeros(nIter,1);
            g_global_size_sumsquared = zeros(nIter,1);
            g_global_fromtruth_sum = zeros(nIter,1);
            g_global_fromtruth_sumsquared = zeros(nIter,1);
        end
        
        
        g_global_count_sum = g_global_count_sum + global_count_sum;
        g_global_count_sumsquared = g_global_count_sumsquared + global_count_sumsquared;
        g_global_coverage_sum = g_global_coverage_sum + global_coverage_sum;
        g_global_coverage_sumsquared = g_global_coverage_sumsquared + global_coverage_sumsquared;
        g_global_coverage_cum_sum = g_global_coverage_cum_sum + global_coverage_cum_sum;
        g_global_coverage_cum_sumsquared = g_global_coverage_cum_sumsquared + global_coverage_cum_sumsquared;
        g_global_speed_sum = g_global_speed_sum + global_speed_sum;
        g_global_speed_sumsquared = g_global_speed_sumsquared + global_speed_sumsquared;
        g_global_move_sum = g_global_move_sum + global_move_sum;
        g_global_move_sumsquared = g_global_move_sumsquared + global_move_sumsquared;
        g_global_size_sum = g_global_size_sum + global_size_sum;
        g_global_size_sumsquared = g_global_size_sumsquared + global_size_sumsquared;
        g_global_fromtruth_sum = g_global_fromtruth_sum + global_fromtruth_sum;
        g_global_fromtruth_sumsquared = g_global_fromtruth_sumsquared + global_fromtruth_sumsquared;     

    end
    
    CI_INT = 0.95/2 + 0.5;
    N = validFiles;
    df = N - 1;

    % Computing g_global stats
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
    
    dataFileName = [DUMPDIR 'clusters_macro_avg.csv'];
    write_csv_headers(dataFileName, headers_clusters_macro_avg);
    fidClustersMacroAvg = fopen(dataFileName,'a');
    
    % Saving all iterations
    for z = 1:nIter
        clu_macro_avg_string = sprintf('"%s",%u,%.4,f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f', ...
            subDir, totalFiles, t_count_avg(z), t_count_sd(z), t_count_se(z), t_count_ci(z), ...
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



