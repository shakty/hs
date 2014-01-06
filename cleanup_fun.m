function cleanup_fun(path2sim)

    % Deletes all temporary computations making by ANAL and MERGE.
    % Takes a path to a root folder, and then goes in each directory
    % (besides aggr and img) and deletes dirs: clusters, agents,
    % truthradius.

    tic
    
    % Every subdirectory of path2sim contains simulations results.
    dirs = dir(path2sim);
    dirIndex = find([dirs.isdir]);

    if (isempty(dirIndex))
        error('Invalid Directory Selected');
    end

    dirs = dirs(dirIndex);
    
    matlabpool(2)
    
    % Each subdir containing partials results will be deleted
    parfor d = 1:length(dirs)

        subDir = dirs(d).name;

        if ( ...
            strcmpi(subDir,'.') ...
            || strcmpi(subDir,'..') ...
            || strcmpi(subDir,'aggr') ...
            || strcmpi(subDir,'img') ...        
        )
            continue;
        end

        dirPath = [path2sim subDir '/'];

        if (exist([dirPath 'agents'], 'dir') == 1)
            rmdir([dirPath 'agents'], 's');
        end
        
        if (exist([dirPath 'clusters'], 'dir') == 1)
            rmdir([dirPath 'clusters'], 's');
        end
        
        if (exist([dirPath 'truthradius'], 'dir') == 1)
            rmdir([dirPath 'truthradius'], 's');
        end

    end
    
    matlabpool close
    
    toc;
end