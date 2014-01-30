function LSF_analysis(path2conf)

% function LSF_analysis( DUMPDIR, simName, sigmas, DUMP, DUMP_RATE, ...
%     RADIUSs, STAY_FOR, CONSENSUS_ON_TRUTH_FOR, CONSENSUS_THRESHOLD, ... % Radius
%     PRECISION, ... % Agents
%     CLU_CUTOFF ... % Clusters
% )
    tic
    
    % Loads all variables that are described in the commented method.
    load([path2conf 'params_all']);
    
    FILES4TASK = 10;
    
    path(path,'util/'); % Help functions
    path(path,'lib/'); % Help functions

    % Import LSF Profile
    parallel.importProfile('/cluster/apps/matlab/support/BrutusLSF8h.settings');
    sched = findResource('scheduler','type','lsf');
    
    % In a dump directory simulations are divided by level of noise (sigma)
    % In each folder (sigma) we start the truthradius analysis.
    for s = 1 : length(sigmas)
        
        sigma = num2str(sigmas(s)*100);
        sigmaSimName = [simName '_s' sigma '/'];
        dumpDir = [DUMPDIR simName '/' sigmaSimName];
        
        files = dir(dumpDir);
        fileIndex = find(~[files.isdir]);

        if (isempty(fileIndex))
            error(['Invalid Directory Selected: ' dumpDir]);
        end
        
        % Log: Matlab will try to create intermediate non-existing folders.
        logFolder = ['log/' simName];
        if (exist(logFolder, 'dir') == 0)
            mkdir(logFolder);
        end
        
        % Dump Folder        
        submitArgs = [' -W 8:00 -R "rusage[mem=8000]" -o ' logFolder '/' simName '.log'];
        set(sched, 'SubmitArguments', submitArgs);
        set(sched, 'DataLocation', [logFolder '/']);
        
        % Create a new Job for each sub folder.
        j = createJob(sched, 'name', ['analysis_' sigma]);
        
        paramsArgs = cell(FILES4TASK, 1);

        outDirAgents = [dumpDir '/' 'agents/'];
        outDirRadius = [dumpDir '/' 'truthradius/'];
        outDirClusters = [dumpDir '/' 'clusters/'];
        
        % Creating outDir if not existing.
        if (exist(outDirAgents, 'dir') == 0)
            mkdir(outDirAgents);
        end
        
        % Creating outDir if not existing.
        if (exist(outDirRadius, 'dir') == 0)
            mkdir(outDirRadius);
        end
        
        % Creating outDir if not existing.
        if (exist(outDirClusters, 'dir') == 0)
            mkdir(outDirClusters);
        end
        
        % Number of files.
        nFiles = length(fileIndex);
        validFiles = 0;
        
        % Load all parameters matrices in one.
        for f = 1:nFiles

            append = files(fileIndex(f)).name;
            fileName = [dumpDir, append];

            % We load only .mat
            [PATH, NAME, EXT] = fileparts(fileName);            
            if (~strcmpi(EXT,'.mat') || ~isempty(strfind(NAME, 'sums_')))  
                continue;
            end
            
            validFiles = validFiles + 1;
            
            % Extracting the part 1-1 from 1-1.mat
            % 6 = length('sums_') + 1; 4 = length('.mat');
            simnameidx = NAME(1:length(NAME));

            paramsObj = struct( ...
                    'folderName', [DUMPDIR simName '/'], ...
                    'simName', sigmaSimName, ...
                    'fileName', simnameidx, ...
                    'RADIUSs', RADIUSs, ...
                    'STAY_FOR', STAY_FOR, ...
                    'CONSENSUS_ON_TRUTH_FOR', CONSENSUS_ON_TRUTH_FOR, ...
                    'CONSENSUS_THRESHOLD', CONSENSUS_THRESHOLD, ...
                    'DUMP', DUMP_ANALYSIS, ...
                    'DUMP_RATE', DUMP_RATE_ANALYSIS, ...
                    'PLOTS', 0, ...
                    'CLU_CUTOFF', CLU_CUTOFF, ...
                    'PRECISION', PRECISION, ...
                    'outDirRadius', outDirRadius, ...
                    'outDirAgents', outDirAgents, ...
                    'outDirClusters', outDirClusters ...      
            );
            
            idx = mod(validFiles, FILES4TASK);
            
            % paramsArgs needs to be enclosed in two cells because
            % otherwise matlab thinks that he needs to pass each cell as a
            % input parameter of wrapperanalysis instead of giving it all
            % to the function that will loop through them.            
            if (idx == 0)
                paramsArgs{FILES4TASK} = paramsObj;
                createTask(j, @wrapperanalysis, 0, {{paramsArgs}});
                paramsArgs = cell(FILES4TASK, 1);
            else
                paramsArgs{idx} = paramsObj;
            end
            
        end % File loop

        % Add files that are left over.
        if (idx ~= 0)
            createTask(j, @wrapperanalysis, 0, {{paramsArgs}});
        end
        
        % Diary works only with batch.
        % diary(j, [logFolder 'job_' sigma '.diary']);
        % Submit all the truthradius tasks.
        submit(j);

        % Create a task that waits to aggregate their results.
%         paramsObjWait = paramsObj;
%         paramsObjWait.j = j;
%         paramsObj.LSF = 1;
%         jWait = createJob(sched, 'name', 'wait2analize');
%         createTask(jWait, @aggregate_onesim, 0, {paramsObjWait});
%         submit(jWait);

    end % Sigmas loop
    toc 
end

%%




