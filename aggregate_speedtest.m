function aggregate_speedtest(dumpDir, subDir, outDir)
    tic;
    
    % Creating outDir if not existing.
    if (exist(outDir, 'dir') == 0)
        mkdir(outDir);
    end
  
    
    % AGENTS
    
    dirPath = [dumpDir subDir ];
    
    files = dir(dirPath);
    fileIndex = find(~[files.isdir]);

    if (isempty(fileIndex))
        error('Invalid Directory Selected');
    end

    % Number of files
    nFiles = length(fileIndex);    
 
    speedTestFileName = [outDir 'speedtest.csv'];
    
    fidSpeedTest = fopen(speedTestFileName,'a'); 
    
    % Scan all files
    for f = 1:nFiles

        fileName = files(fileIndex(f)).name;
        path2file = [dirPath fileName];

        % We load only sums_**.mat (not sums_all)
        [PATH, NAME, EXT] = fileparts(path2file);
        if (~strcmpi(EXT, '.mat'))
            continue;
        end
       
        load(path2file);
        
        p = dump.parameters;
        row = sprintf('"%s",%i,%u,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%i,%i,%.4f,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i', ... 
            subDir, ...
            dump.sim, ...
            dump.run, ...        
            p.R, ...
            p.alpha, ...
            p.tau, ...
            p.v_scaling, ...
            p.sigma, ...
            p.epsilon, ...
            p.seed, ...
            p.nof_cluster, ...
            p.clustersInCircleOfRadius, ...
            dump.consensus1, ...
            dump.consensus10, ...
            dump.consensus20, ...              
            dump.consensus25, ...
            dump.consensus40, ...
            dump.consensus50, ...              
            dump.consensus60, ...
            dump.consensus70, ...
            dump.consensus75, ...           
            dump.consensus80, ...
            dump.consensus90, ...
            dump.consensus100, ...
            dump.ccount1, ...
            dump.ccount10, ...
            dump.ccount20, ...              
            dump.ccount25, ...
            dump.ccount40, ...
            dump.ccount50, ...              
            dump.ccount60, ...
            dump.ccount70, ...
            dump.ccount75, ...           
            dump.ccount80, ...
            dump.ccount90, ...
            dump.ccount100 ...
        );
        fprintf(fidSpeedTest, '%s\n', row);
        
       
        clearvars dump;
        
    end
   
    fclose(fidSpeedTest);

end
