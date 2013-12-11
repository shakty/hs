%% Aggregates the results of the analysis of the simulation results.
tic;

close all;
clear;
clc;

%% Add other directories to path
path(path,'util/'); % Help functions
path(path,'lib/'); % Help functions

aggrParams = 1;
RADIUSs = [.01, .05, .1, .25, .4];
nRadiusesPlusOne = length(RADIUSs) + 1;

DUMPDIR = '/home/stefano/hs/test/';

simName = 'NEWTEST-2013-12-8-17-49/';

path2sim = [DUMPDIR simName];

% Every subdirectory of path2sim contains simulations results.
dirs = dir(path2sim);
dirIndex = find([dirs.isdir]);

if (isempty(dirIndex))
    error('Invalid Directory Selected');
end

outDir = [path2sim 'aggr/'];

% Creating outDir if not existing.
if (exist(outDir, 'dir') == 0)
    mkdir(outDir);
end

agentsFileName = [outDir 'agents.csv'];
paramsFileName = [outDir 'params.csv'];
clustersMacroFileName = [outDir 'clusters_macro.csv'];
truthradiusFileName = [outDir 'truthradius.csv'];

validFiles = 0;

% Each subdir containing results must be aggregated (they are divided by
% run) and then the aggregated results must be aggregated overall.
for d = 1:length(dirIndex)

    subDir = dirs(dirIndex(d)).name;

    if ( ...
        strcmpi(subDir,'.') ...
        || strcmpi(subDir,'..') ...
        || strcmpi(subDir,'aggr') ...
        || strcmpi(subDir,'img') ...        
    )
        continue;
    end

    subDir = [subDir '/'];
    
    validFiles = validFiles + 1;

    dirPath = [path2sim subDir];

    % Aggregate the results of each sub-simulation (level of sigma).
    % The aggregate output is saved in the same dir under a folder called
    % agents/, clusters/, truthradius/.
    aggregate_agents(path2sim, subDir, dirPath, aggrParams); %params.
    aggregate_clusters(path2sim, subDir, dirPath);        
    aggregate_truthradius(path2sim, subDir, dirPath, RADIUSs);
    
    % Paths to CSV files.
    agentsFileCSV = [dirPath 'agents/agents.csv'];
    paramsFileCSV = [dirPath 'agents/params.csv'];
    clustersMacroFileCSV = [dirPath 'clusters/clusters_macro.csv'];
    truthradiusFileCSV = [dirPath 'truthradius/truthradius.csv'];    
    
    % Aggregate all levels of sigmas.
    
    % Clusters.
    load([ dirPath 'clusters/sums_all']);
    
    % Init the arrays on first iteration. Also take full csv files with
    % headers. Then the first line will be removed while merging.
    if (validFiles == 1)
        
        nIter = length(g_global_count_sum);
        
        % Agents.
        ag_global_coverage_sum = zeros(nIter,1);
        ag_global_coverage_sumsquared = zeros(nIter,1);
        ag_global_coverage_cum_sum = zeros(nIter,1);
        ag_global_coverage_cum_sumsquared = zeros(nIter,1);
        ag_global_speed_sum = zeros(nIter,1);
        ag_global_speed_sumsquared = zeros(nIter,1);
        ag_global_movs_sum = zeros(nIter,1);
        ag_global_movs_sumsquared = zeros(nIter,1); 
        ag_global_fromtruth_sum = zeros(nIter,1);
        ag_global_fromtruth_sumsquared = zeros(nIter,1);
        ag_global_pdist_sum =  zeros(nIter,1);
        ag_global_pdist_sumsquared =  zeros(nIter,1);

        % Clusters.
        cg_global_count_sum = zeros(nIter,1);
        cg_global_count_sumsquared = zeros(nIter,1);
        cg_global_maxsize_sum = zeros(nIter,1);
        cg_global_maxsize_sumsquared = zeros(nIter,1);
        cg_global_meansize_sum = zeros(nIter,1);
        cg_global_meansize_sumsquared = zeros(nIter,1); 
        cg_global_speed_sum = zeros(nIter,1);
        cg_global_speed_sumsquared = zeros(nIter,1);
        cg_global_movs_sum = zeros(nIter,1);
        cg_global_movs_sumsquared = zeros(nIter,1);
        cg_global_fromtruth_sum = zeros(nIter,1);
        cg_global_fromtruth_sumsquared = zeros(nIter,1);
        cg_global_bigcpdist_sum = zeros(nIter,1);
        cg_global_bigcpdist_sumsquared = zeros(nIter,1);

        % TruthRadius.
        tg_globalRadiusCounts_sum = zeros(nIter,nRadiusesPlusOne);
        tg_globalRadiusCounts_squared = zeros(nIter,nRadiusesPlusOne);
        tg_globalConsensusOnTruth_sum = zeros(nIter,1);
        tg_globalConsensusOnTruth_squared = zeros(nIter,1);
        
        % Merging CSV files.
        mergeCommand = sprintf('cat %s >> %s', agentsFileCSV, agentsFileName);
        system(mergeCommand);
        
        mergeCommand = sprintf('cat %s >> %s', paramsFileCSV, paramsFileName);
        system(mergeCommand);
        
        mergeCommand = sprintf('cat %s >> %s', clustersMacroFileCSV, clustersMacroFileName);
        system(mergeCommand);
        
        mergeCommand = sprintf('cat %s >> %s', truthradiusFileCSV, truthradiusFileName);
        system(mergeCommand);
    else
        
        % Merging CSV files without headers.
        mergeCommand = sprintf('sed -e ''1d'' %s >> %s', agentsFileCSV, agentsFileName);
        system(mergeCommand);
        
        mergeCommand = sprintf('sed -e ''1d'' %s >> %s', paramsFileCSV, paramsFileName);
        system(mergeCommand);
        
        mergeCommand = sprintf('sed -e ''1d'' %s >> %s', clustersMacroFileCSV, clustersMacroFileName);
        system(mergeCommand);
        
        mergeCommand = sprintf('sed -e ''1d'' %s >> %s', truthradiusFileCSV, truthradiusFileName);
        system(mergeCommand);
    end
    
    % Updates stats arrays.
    
    % Clusters (already loaded above).
    cg_global_count_sum = cg_global_count_sum + g_global_count_sum;
    cg_global_count_sumsquared = cg_global_count_sumsquared + g_global_count_sumsquared;
    cg_global_meansize_sum = cg_global_meansize_sum + g_global_meansize_sum;
    cg_global_meansize_sumsquared = cg_global_meansize_sumsquared + g_global_meansize_sumsquared;
    cg_global_maxsize_sum = cg_global_maxsize_sum + g_global_maxsize_sum;
    cg_global_maxsize_sumsquared = cg_global_maxsize_sumsquared + g_global_maxsize_sumsquared;        
    cg_global_speed_sum = cg_global_speed_sum + g_global_speed_sum;
    cg_global_speed_sumsquared = cg_global_speed_sumsquared + g_global_speed_sumsquared;
    cg_global_movs_sum = cg_global_movs_sum + g_global_movs_sum;
    cg_global_movs_sumsquared = cg_global_movs_sumsquared + g_global_movs_sumsquared;
    cg_global_fromtruth_sum = cg_global_fromtruth_sum + g_global_fromtruth_sum;
    cg_global_fromtruth_sumsquared = cg_global_fromtruth_sumsquared + g_global_fromtruth_sumsquared;
    cg_global_bigcpdist_sum = cg_global_bigcpdist_sum + g_global_bigcpdist_sum;
    cg_global_bigcpdist_sumsquared = cg_global_bigcpdist_sumsquared + g_global_bigcpdist_sumsquared;
        
    % Delete file variables to be sure we do not overwrite some by
    % mistake.
    clearvars g_global_*;
    
    % Agents
    load([ dirPath 'agents/sums_all']);
    
    ag_global_coverage_sum = ag_global_coverage_sum + g_global_coverage_sum;
    ag_global_coverage_sumsquared = ag_global_coverage_sumsquared + g_global_coverage_sumsquared;
    ag_global_coverage_cum_sum = ag_global_coverage_cum_sum + g_global_coverage_cum_sum;
    ag_global_coverage_cum_sumsquared = ag_global_coverage_cum_sumsquared + g_global_coverage_cum_sumsquared;
    ag_global_speed_sum = ag_global_speed_sum + g_global_speed_sum;
    ag_global_speed_sumsquared = ag_global_speed_sumsquared + g_global_speed_sumsquared;
    ag_global_movs_sum = ag_global_movs_sum + g_global_movs_sum;
    ag_global_movs_sumsquared = ag_global_movs_sumsquared + g_global_movs_sumsquared;
    ag_global_fromtruth_sum = ag_global_fromtruth_sum + g_global_fromtruth_sum;
    ag_global_fromtruth_sumsquared = ag_global_fromtruth_sumsquared + g_global_fromtruth_sumsquared;     
    ag_global_pdist_sum = ag_global_pdist_sum + g_global_pdist_sum;
    ag_global_pdist_sumsquared = ag_global_pdist_sumsquared + g_global_pdist_sumsquared;

    % Delete file variables to be sure we do not overwrite some by
    % mistake.
    clearvars g_global_*;
     
    % TruthRadius.
    load([ dirPath 'truthradius/sums_all']);
    
    tg_globalRadiusCounts_sum = tg_globalRadiusCounts_sum + g_globalRadiusCounts_sum;
    tg_globalRadiusCounts_squared = tg_globalRadiusCounts_squared + g_globalRadiusCounts_squared;
    tg_globalConsensusOnTruth_sum = tg_globalConsensusOnTruth_sum + g_globalConsensusOnTruth_sum;
    tg_globalConsensusOnTruth_squared = tg_globalConsensusOnTruth_squared + g_globalConsensusOnTruth_squared;
    
    clearvars g_global*;
    
end


% GLOBAL GLOBAL STATS

    % N: observations.
    N = validFiles;
    % df: degree of freedom.
    df = N - 1;
    % 95 percent confidence interval.
    CI_INT = 0.95/2 + 0.5;
    
    % Clusters
    
    t_count_avg = cg_global_count_sum / N;
    t_count_sd = sqrt(((cg_global_count_sumsquared - ((cg_global_count_sum).^2 / N))) / df);
    t_count_se = t_count_sd / sqrt(N);
    t_count_ci = t_count_se * tquant(CI_INT, df);
    
    t_meansize_avg = cg_global_meansize_sum / N; 
    t_meansize_sd = sqrt(((cg_global_meansize_sumsquared - ((cg_global_meansize_sum).^2 / N))) / df);
    t_meansize_se = t_meansize_sd / sqrt(N);  
    t_meansize_ci = t_meansize_se * tquant(CI_INT, df);
    
    t_maxsize_avg = cg_global_maxsize_sum / N; 
    t_maxsize_sd = sqrt(((cg_global_maxsize_sumsquared - ((cg_global_maxsize_sum).^2 / N))) / df);
    t_maxsize_se = t_meansize_sd / sqrt(N);  
    t_maxsize_ci = t_meansize_se * tquant(CI_INT, df);
    
    t_speed_avg = cg_global_speed_sum / N; 
    t_speed_sd = sqrt(((cg_global_speed_sumsquared - ((cg_global_speed_sum).^2 / N))) / df);
    t_speed_se = t_speed_sd / sqrt(N);  
    t_speed_ci = t_speed_se * tquant(CI_INT, df);
    
    t_move_avg = cg_global_movs_sum / N; 
    t_move_sd = sqrt(((cg_global_movs_sumsquared - ((cg_global_movs_sum).^2 / N))) / df);
    t_move_se = t_move_sd / sqrt(N);  
    t_move_ci = t_move_se * tquant(CI_INT, df);

    t_fromtruth_avg = cg_global_fromtruth_sum / N; 
    t_fromtruth_sd = sqrt(((cg_global_fromtruth_sumsquared - ((cg_global_fromtruth_sum).^2 / N))) / df);
    t_fromtruth_se = t_fromtruth_sd / sqrt(N);  
    t_fromtruth_ci = t_fromtruth_se * tquant(CI_INT, df);
    
    t_bigcpdist_avg = cg_global_bigcpdist_sum / N; 
    t_bigcpdist_sd = sqrt(((cg_global_bigcpdist_sumsquared - ((cg_global_bigcpdist_sum).^2 / N))) / df);
    t_bigcpdist_se = t_bigcpdist_sd / sqrt(N);
    t_bigcpdist_ci = t_bigcpdist_se * tquant(CI_INT, df);
    
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
        clu_macro_avg_string = sprintf('"%s",%u,%u,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f', ...
            simName, N, z, ...
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
    
    % Agents    
    
    t_cover_avg = ag_global_coverage_sum / N; 
    t_cover_sd = sqrt(((ag_global_coverage_sumsquared - ((ag_global_coverage_sum).^2 / N))) / df);
    t_cover_se = t_cover_sd / sqrt(N);  
    t_cover_ci = t_cover_se * tquant(CI_INT, df);
    
    t_cover_cum_avg = ag_global_coverage_cum_sum / N; 
    t_cover_cum_sd = sqrt(((ag_global_coverage_cum_sumsquared - ((ag_global_coverage_cum_sum).^2 / N))) / df);
    t_cover_cum_se = t_cover_cum_sd / sqrt(N);  
    t_cover_cum_ci = t_cover_cum_se * tquant(CI_INT, df);
    
    t_speed_avg = ag_global_speed_sum / N; 
    t_speed_sd = sqrt(((ag_global_speed_sumsquared - ((ag_global_speed_sum).^2 / N))) / df);
    t_speed_se = t_speed_sd / sqrt(N);  
    t_speed_ci = t_speed_se * tquant(CI_INT, df);
    
    t_move_avg = ag_global_movs_sum / N; 
    t_move_sd = sqrt(((ag_global_movs_sumsquared - ((ag_global_movs_sum).^2 / N))) / df);
    t_move_se = t_move_sd / sqrt(N);  
    t_move_ci = t_move_se * tquant(CI_INT, df);
 
    t_fromtruth_avg = ag_global_fromtruth_sum / N; 
    t_fromtruth_sd = sqrt(((ag_global_fromtruth_sumsquared - ((ag_global_fromtruth_sum).^2 / N))) / df);
    t_fromtruth_se = t_fromtruth_sd / sqrt(N);
    t_fromtruth_ci = t_fromtruth_se * tquant(CI_INT, df);
    
    t_pdist_avg = ag_global_pdist_sum / N; 
    t_pdist_sd = sqrt(((ag_global_pdist_sumsquared - ((ag_global_pdist_sum).^2 / N))) / df);
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
            simName, N, z, ...
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
    
    % Truthradius

        % Computing g_global stats.
    
     headers_truthradius_avg = {
        'simname', ...
        'N', ...
        't' ...
    };
    
    avgStartFrom = length(headers_truthradius_avg) + 1;
    
    row_string_data_cell = cell(nRadiusesPlusOne*4, 1);
    
    % Computing global stats            
    for i = 1 : nRadiusesPlusOne
        
        % This part is actally different from aggregate_truthradius, that
        % has cell array, and r_out is in the radius array.
        % Not a big deal anyway.
        if (i < nRadiusesPlusOne)            
            radiusStr = ['r_' strrep(num2str(RADIUSs(i)), '0.', '')];
        else
            radiusStr = 'r_out';
        end
        
        hIdx = (i-1)*4 + avgStartFrom;
        meanName = ['t_' radiusStr '_mean'];      
        eval([meanName ' = tg_globalRadiusCounts_sum(:,i) / N;']);

        sdName = ['t_' radiusStr '_sd'];        
        % sdValue = sqrt(((g_globalRadiusCounts_squared(i) - ((g_globalRadiusCounts_sum(i)).^2 / N))) / df);
        eval([sdName ' = sqrt(((tg_globalRadiusCounts_squared(:,i) - ((tg_globalRadiusCounts_sum(:,i)).^2 / N))) / df);']);

        seName = ['t_' radiusStr '_se'];
        eval([seName ' = ' sdName ' / sqrt(N);']);

        ciName = genvarname(['t_' radiusStr '_ci']);            
        eval([ciName ' = ' seName ' * tquant(CI_INT, df);']);    
        
        headers_truthradius_avg{hIdx} = regexprep(meanName, 't_', '', 'once');
        headers_truthradius_avg{hIdx + 1} = regexprep(sdName, 't_', '', 'once');
        headers_truthradius_avg{hIdx + 2} = regexprep(seName, 't_', '', 'once');
        headers_truthradius_avg{hIdx + 3} = regexprep(ciName, 't_', '', 'once'); 
        
        hIdx = hIdx - avgStartFrom + 1;
        row_string_data_cell{hIdx} = [meanName '(z)'];
        row_string_data_cell{hIdx + 1} = [sdName '(z)'];
        row_string_data_cell{hIdx + 2} = [seName '(z)'];
        row_string_data_cell{hIdx + 3} = [ciName '(z)'];                
    end
    
    t_consensusOnTruth_avg = tg_globalConsensusOnTruth_sum / N; 
    t_consensusOnTruth_sd = sqrt(((tg_globalConsensusOnTruth_squared - ((tg_globalConsensusOnTruth_sum).^2 / N))) / df);
    t_consensusOnTruth_se = t_consensusOnTruth_sd / sqrt(N);  
    t_consensusOnTruth_ci = t_consensusOnTruth_se * tquant(CI_INT, df);
    
    hIdx = i*4 + avgStartFrom;
    headers_truthradius_avg{hIdx} = 'consensus.avg';
    headers_truthradius_avg{hIdx + 1} = 'consensus.sd';
    headers_truthradius_avg{hIdx + 2} = 'consensus.se';
    headers_truthradius_avg{hIdx + 3} = 'consensus.ci';

    truthradiusAvgFileName = [outDir 'truthradius_avg_all.csv'];
    write_csv_headers(truthradiusAvgFileName, headers_truthradius_avg);
    fidTruthRadiusAvg = fopen(truthradiusAvgFileName, 'a');

    % Saving all iterations
    for z = 1:nIter
        truthradius_avg_string = sprintf('"%s",%u,%u,', simName, N, z);
        for j = 1 : length(row_string_data_cell)
            truthradius_avg_string = [truthradius_avg_string sprintf('%.4f',(eval(row_string_data_cell{j}))) ','];
        end
        truthradius_avg_string = sprintf('%s%.4f,%.4f,%.4f,%.4f', ...
            truthradius_avg_string, ...
            t_consensusOnTruth_avg(z), t_consensusOnTruth_sd(z), ...
            t_consensusOnTruth_se(z), t_consensusOnTruth_ci(z) ...
        );
        fprintf(fidTruthRadiusAvg, '%s\n', truthradius_avg_string);
    end

    fclose(fidTruthRadiusAvg);
    
    toc;