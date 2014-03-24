function speedtest_fun(path2sim)

    tic
    
    path(path,'util/'); % Help functions
    path(path,'lib/'); % Help functions

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
        'R', ...
        'alpha', ...
        'tau', ...
        'init.vscaling', ...
        'sigma', ...
        'epsilon', ...
        'seed', ...
        'init.ccount', ...
        'init.placement', ...
        'init.band.i', ...
        'init.band.y', ...
        'consensus1', ...
        'consensus10', ...
        'consensus20', ...              
        'consensus25', ...
        'consensus30', ...
        'consensus40', ...
        'consensus50', ...              
        'consensus60', ...
        'consensus70', ...
        'consensus75', ...           
        'consensus80', ...
        'consensus90', ...
        'consensus100', ...
        'ccount1', ...
        'ccount10', ...
        'ccount20', ...              
        'ccount25', ...
        'ccount30', ...
        'ccount40', ...
        'ccount50', ...              
        'ccount60', ...
        'ccount70', ...
        'ccount75', ...           
        'ccount80', ...
        'ccount90', ...
        'ccount100', ...
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