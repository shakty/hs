function aggregate_truthradius(params)

    j = params.j;
    outDir = params.outDir;
    
    CI_INT = 0.95/2 + 0.5;

    % Creating variables used to save results to file.
    
    headers_truthradius_avg = {
        'simname', ...
        'simcount', ...
        't' ...
    };

    avg_symbols = { '_mean', '_sd', '_se', '_ci' };
    nSyms = length(avg_symbols);

    % nRadiuses
    nRadiuses = length(RADIUSs);
    hAvgStartFrom = 3;
    radiusesStr = {};
    for i = 1 : nRadiuses
        radiusStr = ['r_' num2str(RADIUSs(i))];
        radiusesStr{i} = radiusStr;
        for j = 1 : nSyms
            idxJ = (i -1) * nSyms + j + hAvgStartFrom;
            radiusesStrAvg = [radiusStr avg_symbols{j}];
            headers_truthradius_avg{idxJ} = radiusesStrAvg;
        end
        
    end
    
    % Adding not in radius
    i = i + 1;
    radiusStr = 'r_out';
    radiusesStr{i} = radiusStr;
    
    for j = 1 : nSyms
        idxJ = (i -1) * nSyms + j + hAvgStartFrom;
        headers_truthradius_avg{idxJ} = [radiusStr avg_symbols{j}];
    end
    
    % Adding flag consensus on Truth.
    i = i + 1;
    radiusStr = 'consensus';
    radiusesStr{i} = radiusStr;
    
    for j = 1 : nSyms
        idxJ = (i -1) * nSyms + j + hAvgStartFrom;
        headers_truthradius_avg{idxJ} = [radiusStr avg_symbols{j}];
    end
    
    % +1 is NOT_IN_RADIUS.
    nRadiusesPlusOne = nRadiuses + 1;
    
    % End creating variables used to save results to file.
    
    
    % GLOBAL variables
  
    globalRadiusCounts = zeros(nIter, nRadiusesPlusOne);
    globalRadiusCounts_squared = zeros(nIter, nRadiusesPlusOne);
    globalConsensusOnTruth = zeros(1, nIter);
    globalConsensusOnTruth_squared = zeros(1, nIter);
    
    wait(j);
    
    % Delete Job
    delete(j);
    
    
    files = dir(outDir);
    fileIndex = find(~[files.isdir]);

    if (isempty(fileIndex))
        error('Invalid Directory Selected');
    end

    validFileIdx = 0;
    summaryObj = {};
    
    % Retrieving global values for the set of simulations from the first one
    
    fileName = [ourDir '1-1.mat'];    
    load(fileName);
    v = dump.agentsv;
    % the total number of time steps per run (assumed constant) and dump rate
    nIter = size(v,3);
    % idxsIters = find(mod(1:nIter, DUMP_RATE) == 0);
    % number of files
    nFiles = length(fileIndex);
    
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
    end
    
    
    
    N = validFileIdx; % last value
    df = N - 1;
    
    
    
    
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
            clu_macro_avg_string = sprintf('"%s",%u,%u,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f', ...
                simName, N, z, ...
                t_cover_avg(z), t_cover_sd(z), t_cover_se(z), t_cover_ci(z), ...
                t_cover_cum_avg(z), t_cover_cum_sd(z), t_cover_cum_se(z), t_cover_cum_ci(z), ...
                t_speed_avg(z), t_speed_sd(z), t_speed_se(z), t_speed_ci(z), ...
                t_movs_avg(z), t_movs_sd(z), t_movs_se(z), t_movs_ci(z), ...
                t_fromtruth_avg(z), t_fromtruth_sd(z), t_fromtruth_se(z), t_fromtruth_ci(z));
            
            fprintf(fidTruthRadiusAvg,'%s\n', clu_macro_avg_string);   
         end

        fclose(fidTruthRadiusAvg);
    end
    
    % Delete the directory outDir

end

