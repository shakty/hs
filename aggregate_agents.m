function aggregate_agents(dumpDir, subDir, outDir, aggrParams)
    tic;
    
    % This function aggregates the params as well.
    
    % Creating outDir if not existing.
    if (exist(outDir, 'dir') == 0)
        mkdir(outDir);
    end
    
    % PARAMS (for all)
    
    if (aggrParams)
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
    end
    
    % AGENTS
    
    dirPath = [dumpDir subDir 'agents/'];
    
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
        agentsFileCSV = [dirPath 'agents_' simnameidx '.csv'];
        mergeCommand = sprintf('cat %s >> %s', ...
            agentsFileCSV, agentsFileName);        
        system(mergeCommand);
        
        if (aggrParams)
            % Just in this case, merge also the params file.
            paramsFileCSV = [dirPath 'params_' simnameidx '.csv'];
            mergeCommand = sprintf('cat %s >> %s', ...
                paramsFileCSV, paramsFileName);        
            system(mergeCommand,'-echo');
        end
        
        % Load file to compute round statistics.
        load(fileName);
        
        % Increment index of valid files.
        validFileIdx = validFileIdx + 1;
        
        % Updates stats arrays.
        g_global_coverage_sum = g_global_coverage_sum + global_coverage_sum;
        g_global_coverage_sumsquared = g_global_coverage_sumsquared + global_coverage_sumsquared;
        g_global_coverage_cum_sum = g_global_coverage_cum_sum + global_coverage_cum_sum;
        g_global_coverage_cum_sumsquared = g_global_coverage_cum_sumsquared + global_coverage_cum_sumsquared;
        g_global_speed_sum = g_global_speed_sum + global_speed_sum;
        g_global_speed_sumsquared = g_global_speed_sumsquared + global_speed_sumsquared;
        g_global_movs_sum = g_global_movs_sum + global_movs_sum;
        g_global_movs_sumsquared = g_global_movs_sumsquared + global_movs_sumsquared;
        g_global_fromtruth_sum = g_global_fromtruth_sum + global_fromtruth_sum;
        g_global_fromtruth_sumsquared = g_global_fromtruth_sumsquared + global_fromtruth_sumsquared;     
        g_global_pdist_sum = g_global_pdist_sum + global_pdist_sum;
        g_global_pdist_sumsquared = g_global_pdist_sumsquared + global_pdist_sumsquared;
        
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
    
    t_move_avg = g_global_movs_sum / N; 
    t_move_sd = sqrt(((g_global_movs_sumsquared - ((g_global_movs_sum).^2 / N))) / df);
    t_move_se = t_move_sd / sqrt(N);  
    t_move_ci = t_move_se * tquant(CI_INT, df);
 
    t_fromtruth_avg = g_global_fromtruth_sum / N; 
    t_fromtruth_sd = sqrt(((g_global_fromtruth_sumsquared - ((g_global_fromtruth_sum).^2 / N))) / df);
    t_fromtruth_se = t_fromtruth_sd / sqrt(N);
    t_fromtruth_ci = t_fromtruth_se * tquant(CI_INT, df);
    
    t_pdist_avg = g_global_pdist_sum / N; 
    t_pdist_sd = sqrt(((g_global_pdist_sumsquared - ((g_global_pdist_sum).^2 / N))) / df);
    t_pdist_se = t_pdist_sd / sqrt(N);  
    t_pdist_ci = t_pdist_se * tquant(CI_INT, df);
    
    headers_agents_avg = {
        'simname', ...
        'N', ...
        't', ...
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
        'fromtruth.avg', ...
        'fromtruth.sd', ...
        'fromtruth.se', ...
        'fromtruth.ci', ...
        'pdist.avg', ...
        'pdist.sd', ...
        'pdist.se', ...
        'pdist.ci' ...        
    };

    agentsAvgFileName = [outDir 'agents_avg_all.csv'];
    write_csv_headers(agentsAvgFileName, headers_agents_avg);
    fidAgentsAvg = fopen(agentsAvgFileName, 'a');

    % Saving all iterations
    for z = 1:nIter
        agents_avg_string = sprintf('"%s",%u,%u,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f', ...
            subDir, N, z, ...
            t_cover_avg(z), t_cover_sd(z), t_cover_se(z), t_cover_ci(z), ...
            t_cover_cum_avg(z), t_cover_cum_sd(z), t_cover_cum_se(z), t_cover_cum_ci(z), ...
            t_speed_avg(z), t_speed_sd(z), t_speed_se(z), t_speed_ci(z), ...
            t_move_avg(z), t_move_sd(z), t_move_se(z), t_move_ci(z), ...            
            t_fromtruth_avg(z), t_fromtruth_sd(z), t_fromtruth_se(z), t_fromtruth_ci(z), ...
            t_pdist_avg(z), t_pdist_sd(z), t_pdist_se(z), t_pdist_ci(z) ...
        );

        fprintf(fidAgentsAvg, '%s\n', agents_avg_string);   
    end

    fclose(fidAgentsAvg);
    
    % Delete the files ?
end
