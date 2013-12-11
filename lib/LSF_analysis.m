function LSF_analysis( DUMPDIR, simName, sigmas, DUMP, DUMP_RATE, ...
    RADIUSs, STAY_FOR, CONSENSUS_ON_TRUTH_FOR, CONSENSUS_THRESHOLD, ... % Radius
    PRECISION, ... % Agents
    CLU_CUTOFF ... % Clusters
)


    FILES4TASK = 10; 
    FILES4TASK_PLUSONE = FILES4TASK + 1;
    
    path(path,'util/'); % Help functions
    path(path,'lib/'); % Help functions

    % Import LSF Profile
    parallel.importProfile('/cluster/apps/matlab/support/BrutusLSF8h.settings');
    sched = findResource('scheduler','type','lsf');
    
    % In a dump directory simulations are divided by level of noise (sigma)
    % In each folder (sigma) we start the truthradius analysis.
    for s = 1 : length(sigmas)
        
        sigma = sigmas(s);
        dumpDir = [DUMPDIR simName '_' sigma '/'];
        outDir = [dumpDir 'tmp/'];

        files = dir(dumpDir);
        fileIndex = find(~[files.isdir]);

        if (isempty(fileIndex))
            error('Invalid Directory Selected');
        end

        % Creating outDir if not existing.
        if (exist(outDir, 'dir') == 0)
            mkdir(outDir);
        end
    
        % Log: Matlab will try to create intermediate non-existing folders.
        logFolder = ['log/' simName];
        mkdir(logFolder);
        
        % Dump Folder        
        submitArgs = [' -R "rusage[mem=8000]" -o ' logFolder '/' simName '.log'];
        set(sched, 'SubmitArguments',submitArgs);
        set(sched, 'DataLocation', [logFolder '/']);
    
        
        % Create a new Job for each sub folder.
        j = createJob(sched, 'name', ['analysis_' sigma]);
        
        % Number of files.
        nFiles = length(fileIndex);
        
        paramsArgs = cell(FILES4TASK, 1);

        % Load all parameters matrices in one.
        for f = 1:nFiles

            append = files(fileIndex(f)).name;
            fileName = [dumpDir, append];

            % We load only .mat
            [PATH, NAME, EXT_tmp] = fileparts(fileName);
            if (~strcmpi(EXT_tmp,'.mat') ...                
                || strfind(NAME, 'sums_') == 1 ...
                || strcmp(NAME_tmp, 'temporalysis') == 1 ...
                || strcmp(NAME_tmp, 'global_count_sum') == 1 ...
                || strcmp(NAME_tmp, 'global_count_sumsquared') == 1 ...
                || strcmp(NAME_tmp, 'radiusCounts') == 1 ...
                || strcmp(NAME_tmp, 'radiusCountsAvg') == 1) ...

                continue;
            end

            paramsObj = struct( ...
                    'folderName', DUMPDIR, ...
                    'simName', simName, ...
                    'fileName', NAME, ...
                    'RADIUSs', RADIUSs, ...
                    'STAY_FOR', STAY_FOR, ...
                    'CONSENSUS_ON_TRUTH_FOR', CONSENSUS_ON_TRUTH_FOR, ...
                    'CONSENSUS_THRESHOLD', CONSENSUS_THRESHOLD, ...
                    'DUMP', DUMP, ...
                    'DUMP_RATE', DUMP_RATE, ...
                    'PLOTS', 0, ...
                    'CLU_CUTOFF', CLU_CUTOFF, ...
                    'PRECISION', PRECISION, ...
                    'outDirRadius', outDirRadius, ...
                    'outDirAgents', outDirAgents, ...
                    'outDirClusters', outDirClusters ...      
            );
            
            idx = mod(f, FILES4TASK_PLUSONE);
            paramsArgs{idx} = paramsObj;
            
            if (idx == FILES4TASK)
                createTask(j, @wrapperanalysis, 0, {paramsArgs});
                paramsArgs = cell(FILES4TASK, 1);
            end
        end % File loop

        % Add files that are left over.
        if (idx ~= FILES4TASK)
            createTask(j, @wrapperanalysis, 0, {paramsArgs});         
        end
        
        % Submit all the truthradius tasks.
        submit(j);

        % Create a task that waits to aggregate their results.
        paramsObjWait = paramsObj;
        paramsObjWait.j = j;
        paramsObj.LSF = 1;
        jWait = createJob(sched, 'name', 'wait2analize');
        createTask(jWait, @aggregate_onesim, 0, {paramsObjWait});
        submit(jWait);

    end % Sigmas loop
    
end

%%




