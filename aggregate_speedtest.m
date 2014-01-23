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
        row = sprintf('"%s",%u,%u,%u,%.4f,%.4f,%.4f,%.4f,%.4f,%u', ... 
            subDir, ...
            dump.sim, ...
            dump.run, ... 
            dump.counter, ...       
            p.R, ...
            p.alpha, ...
            p.v_scaling, ...
            p.sigma, ...
            p.epsilon, ...
            p.seed ...
        );
        fprintf(fidSpeedTest, '%s\n', row);
        
       
        clearvars dump;
        
    end
   
    fclose(fidSpeedTest);

end
