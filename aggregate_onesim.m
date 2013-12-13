function aggregate_onesim(params)

    %% Aggregates the results of the analysis of one simulation
    
    close all;
    clear;
    clc;

    %% Add other directories to path
    path(path,'util/'); % Help functions
    path(path,'lib/'); % Help functions

    aggrParams = 1;
    
    RADIUSs = params.RADIUSs;

    DUMPDIR = params.DUMPDIR;
    simName = params.simName;

    path2sim = [DUMPDIR simName];
    outDir = [path2sim 'aggr/'];

    % Creating outDir if not existing.
    if (exist(outDir, 'dir') == 0)
        mkdir(outDir);
    end
    
    if (params.LSF)
        c = parcluster;
        j = findJob(c,'name', params.jobName);
        wait(j)
    end
    
    tic;
    
    % Aggregate the results of each sub-simulation (level of sigma).
    % The aggregate output is saved in the same dir under a folder called
    % agents/, clusters/, truthradius/.
    aggregate_agents(path2sim, subDir, dirPath, aggrParams); %params.
    aggregate_clusters(path2sim, subDir, dirPath);        
    aggregate_truthradius(path2sim, subDir, dirPath, RADIUSs);

    % Delete job.
    delete(params.j);
    toc;
end