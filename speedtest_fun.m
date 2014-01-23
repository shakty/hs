function speedtest_fun(path2sim)

    tic
    

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

    
    speedTestFileName = [outDir 'speedtest.csv'];
       
    headers = {
        'simname', ...
        'simcount', ...
        'run', ...
        't', ...
        'R', ...
        'alpha', ...
        'init.vscaling', ...
        'sigma', ...
        'epsilon', ...
        'seed' ...
    };

    
    % This function overwrites exiting files.
    write_csv_headers(speedTestFileName, headers); 
    
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
        
        aggregate_speedtest(path2sim, subDir, outDir)
    end
    
    toc;
end