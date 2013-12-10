function aggregate_agents(dumpDir, subDir, outDir, RADIUSs)
    tic;
    
    % This function aggregates the params as well.
    
    % Creating outDir if not existing.
    if (exist(outDir, 'dir') == 0)
        mkdir(outDir);
    end
    
    % TRUTHRADIUS
    
    % nRadiuses
    nRadiuses = length(RADIUSs);
    nRadiusesPlusOne = nRadiuses + 1;
    
    dirPath = [dumpDir subDir 'truthradius/'];
    
    % PreLoad one file (to create the variable nIter).
    fileName = 'sums_1-1.mat';    
    load([dirPath fileName]);
    nIter = length(globalRadiusCounts);
    
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
    
    g_globalRadiusCounts_sum = zeros(nIter,nRadiusesPlusOne);
    g_globalRadiusCounts_squared = zeros(nIter,nRadiusesPlusOne);
    g_globalConsensusOnTruth_sum = zeros(nIter,1);
    g_globalConsensusOnTruth_squared = zeros(nIter,1); 
    
    % Creating variables used to save results to file.
    headers_truthradius = {
        'simname', ...
        'simcount', ...
        'run', ...
        't' ...
    };

    radiusesStr = cell(nRadiusesPlusOne + 1, 1);
    
    hStartFrom = 4;
    for i = 1 : nRadiuses
        radiusStr = ['r_' strrep(num2str(RADIUSs(i)), '0.', '')];
        radiusesStr{i} = radiusStr;
        headers_truthradius{i + hStartFrom} = radiusStr;
    end
    
    % Adding not in radius
    i = i + 1;
    radiusStr = 'r_out';
    radiusesStr{i} = radiusStr;
    headers_truthradius{i + hStartFrom} = radiusStr;
    
    % Adding flag consensus on Truth.
    i = i + 1;
    radiusStr = 'consensus';
    radiusesStr{i} = radiusStr;
    headers_truthradius{i + hStartFrom} = radiusStr;
   
    truthradiusFileName = [outDir 'truthradius.csv'];
    % This function overwrites exiting files.
    write_csv_headers(truthradiusFileName, headers_truthradius);   
    
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
        truthradiusFileCSV = [dirPath 'truthradius_' simnameidx '.csv'];
        mergeCommand = sprintf('cat %s >> %s', ...
            truthradiusFileCSV, truthradiusFileName);        
        system(mergeCommand);
        
        % Load file to compute round statistics.
        load(fileName);
        
        % Increment index of valid files.
        validFileIdx = validFileIdx + 1;
        
        % Updates stats arrays.
        g_globalRadiusCounts_sum = g_globalRadiusCounts_sum + globalRadiusCounts;
        g_globalRadiusCounts_squared = g_globalRadiusCounts_squared + globalRadiusCounts_squared;
        g_globalConsensusOnTruth_sum = g_globalConsensusOnTruth_sum + globalConsensusOnTruth;
        g_globalConsensusOnTruth_squared = g_globalConsensusOnTruth_squared + globalConsensusOnTruth_squared;
    
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
    
     headers_truthradius_avg = {
        'simname', ...
        'N', ...
        't' ...
    };
    
    avgStartFrom = length(headers_truthradius_avg);
    
    
    row_string_sprintf = '"%s",%u,%u,';
    row_string_data = 'subDir, N, z,';
    row_string_data_cell = cell(nRadiusesPlusOne*4, 1);
    
    % Computing global stats            
    for i = 1 : nRadiusesPlusOne
        hIdx = (i-1)*4 + avgStartFrom;
        meanName = ['t_' radiusesStr{i} '_mean'];        
        eval([meanName ' = g_globalRadiusCounts_sum(i) / N;']);

        sdName = ['t_' radiusesStr{i} '_sd'];        
        % sdValue = sqrt(((g_globalRadiusCounts_squared(i) - ((g_globalRadiusCounts_sum(i)).^2 / N))) / df);
        eval([sdName ' = sqrt(((g_globalRadiusCounts_squared(i) - ((g_globalRadiusCounts_sum(i)).^2 / N))) / df);']);

        seName = ['t_' radiusesStr{i} '_se'];
        eval([seName ' = ' sdName ' / sqrt(N);']);

        ciName = genvarname(['t_' radiusesStr{i} '_ci']);            
        eval([ciName ' = ' seName ' * tquant(CI_INT, df);']);    
        
        headers_truthradius_avg{hIdx} = strrep(meanName, 't_', '');
        headers_truthradius_avg{hIdx + 1} = strrep(sdName, 't_', '');
        headers_truthradius_avg{hIdx + 2} = strrep(seName, 't_', '');
        headers_truthradius_avg{hIdx + 3} = strrep(ciName, 't_', ''); 
        
        hIdx = hIdx - avgStartFrom + 1;
        row_string_data_cell{hIdx} = [meanName '(z)'];
        row_string_data_cell{hIdx + 1} = [sdName '(z)'];
        row_string_data_cell{hIdx + 2} = [seName '(z)'];
        row_string_data_cell{hIdx + 3} = [ciName '(z)'];
        
        row_string_sprintf = [row_string_sprintf ',%u,%u,%u,%u'];
        row_string_data = [row_string_data meanName '(z),' sdName '(z),'];
        row_string_data = [row_string_data seName '(z),' ciName '(z),'];
        
    end
    
    t_consensusOnTruth_avg = g_globalConsensusOnTruth_sum / N; 
    t_consensusOnTruth_sd = sqrt(((g_globalConsensusOnTruth_squared - ((g_globalConsensusOnTruth_sum).^2 / N))) / df);
    t_consensusOnTruth_se = t_consensusOnTruth_sd / sqrt(N);  
    t_consensusOnTruth_ci = t_consensusOnTruth_se * tquant(CI_INT, df);
    
    row_string_sprintf = [row_string_sprintf ',%u,%u,%u,%u'];
    row_string_data = [row_string_data 't_consensusOnTruth_avg(z),t_consensusOnTruth_sd(z),'];
    row_string_data = [row_string_data 't_consensusOnTruth_se(z),t_consensusOnTruth_ci(z)'];
    
            
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
        truthradius_avg_string = sprintf('"%s",%u,%u,', subDir, N, z);
        for j = 1 : length(row_string_data_cell)
            truthradius_avg_string = [truthradius_avg_string sprintf('%.4f',(eval(row_string_data_cell{j}))) ','];
        end
        truthradius_avg_string = sprintf('%s,%.4f,%.4f,%.4f,%.4f', ...
            truthradius_avg_string, t_consensusOnTruth_avg, ...
            t_consensusOnTruth_se, t_consensusOnTruth_ci ...
        );
        % eval([truthradius_avg_string '= [row_string_data_array];']);
        % eval([truthradius_avg_string '= [row_string_data];']);
        % eval([truthradius_avg_string ' = sprintf(row_string_sprintf, row_string_data);']);
        fprintf(fidTruthRadiusAvg, '%s\n', truthradius_avg_string);   
    end

    fclose(fidTruthRadiusAvg);
    
    % Delete the files ?
end
