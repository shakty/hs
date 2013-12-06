function temporal_analysis_agents( DUMPDIR, simName, RADIUSs, STAY_FOR, ...
    CONSENSUS_ON_TRUTH_FOR, CONSENSUS_THRESHOLD, CSV_DUMP, DUMP_RATE, PLOTS)

%     path(path,'util/'); % Help functions
%     path(path,'lib/'); % Help functions
% 
%     % Maybe import it in a higher level
%     parallel.importProfile('/cluster/apps/matlab/support/BrutusLSF8h.settings');
% 
%     % How many tasks group in one job. Should be multiple with SIMS4TASK.
%     TASKS4JOB = 20;
% 
%     jobCount = 1;
% 
%     logFolder = ['log/' params.simName];
%     mkdir(logFolder); % The name is unique under the dump directory.
%     dumpFolder = [ params.dumpDir params.simName];
% 
%     sched = findResource('scheduler','type','lsf');
%     %sched=parcluster('BrutusLSF8h');
%     emails = 'sbalietti@ethz.ch';
%     submitArgs = [' -R "rusage[mem=8000]" -o ' logFolder '/' params.simName '.log'];
%     set(sched, 'SubmitArguments',submitArgs);
%     set(sched, 'DataLocation', [logFolder '/']);
% 
%     j = createJob(sched);
% 
%     jobCount = jobCount + 1;

    jobCell = {};


    % Counts the number of agents that are distant R from the truth
    % R is every value in the array RADIUSs
    
    % To be considered stable the agent must stay in the same within
    % the same radius for STAY_FOR steps
    
    headers_truthradius = {
        'simname', ...
        'simcount', ...
        'run', ...
        't' ...
    };
    
    headers_truthradius_avg = {
        'simname', ...
        'simcount', ...
        't' ...
    };

    avg_symbols = { '_mean', '_sd', '_se', '_ci' };
    nSyms = length(avg_symbols);

    row_string_sprintf = '"%s",%u,%u,';
    
    % nRadiuses
    nRadiuses = length(RADIUSs);
    hStartFrom = 4;
    hAvgStartFrom = 3;
    radiusesStr = {};
    for i = 1 : nRadiuses
        radiusStr = ['r_' num2str(RADIUSs(i))];
        radiusesStr{i} = radiusStr;
        headers_truthradius{i + hStartFrom} = radiusStr;
        row_string_sprintf = [row_string_sprintf ',%u'];
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
    headers_truthradius{i + hStartFrom} = radiusStr;
    row_string_sprintf = [row_string_sprintf ',%u'];
    
    for j = 1 : nSyms
        idxJ = (i -1) * nSyms + j + hAvgStartFrom;
        headers_truthradius_avg{idxJ} = [radiusStr avg_symbols{j}];
    end
    
    % Adding flag consensus on Truth.
    i = i + 1;
    radiusStr = 'consensus';
    radiusesStr{i} = radiusStr;
    headers_truthradius{i + hStartFrom} = radiusStr;
    row_string_sprintf = [row_string_sprintf ',%u']; 
    
    for j = 1 : nSyms
        idxJ = (i -1) * nSyms + j + hAvgStartFrom;
        headers_truthradius_avg{idxJ} = [radiusStr avg_symbols{j}];
    end
    
    % Not in radius.
    NOT_IN_RADIUS = -1;
    % +1 is NOT_IN_RADIUS.
    nRadiusesPlusOne = nRadiuses + 1;
    
    % Truth radius idx.
    TRUTH_RADIUS = 1;
    
     % CI
    CI_INT = 0.95/2 + 0.5;

    dumpDir = [DUMPDIR simName '/'];
    
    files = dir(dumpDir);
    fileIndex = find(~[files.isdir]);

    if (isempty(fileIndex))
        error('Invalid Directory Selected');
    end

    validFileIdx = 0;
    
    % Retrieving global values for the set of simulations from the first one
    
    fileName = [dumpDir '1-1.mat'];    
    load(fileName);
    v = dump.agentsv;
    % the total number of time steps per run (assumed constant) and dump rate
    nIter = size(v,3);
    nAgents = size(v,2);
    AGENTS_4_CONSENSUS = floor(nAgents * CONSENSUS_THRESHOLD);
    
    % number of files
    nFiles = length(fileIndex);
    
    % GLOBAL variables
  
    globalRadiusCounts = zeros(nIter, nRadiusesPlusOne);
    globalRadiusCounts_squared = zeros(nIter, nRadiusesPlusOne);
    globalConsensusOnTruth = zeros(1, nIter);
    globalConsensusOnTruth_squared = zeros(1, nIter);
    
    if (CSV_DUMP)

        dataFileName = [dumpDir 'truth_radius.csv'];
        write_csv_headers(dataFileName, headers_truthradius);
        fidTruthRadius = fopen(dataFileName,'a');
                
        dataFileName = [dumpDir 'truth_radius_avg.csv'];
        write_csv_headers(dataFileName, headers_truthradius_avg);
        fidTruthRadiusAvg = fopen(dataFileName,'a');
       
    end
    
    % Load all parameters matrices in one
    for f = 1:nFiles

        append = files(fileIndex(f)).name;
        fileName = [dumpDir, append];

        % We load only .mat
        [PATH, NAME_tmp, EXT_tmp] = fileparts(fileName);
        if (~strcmpi(EXT_tmp,'.mat') ...
            || strcmp(NAME_tmp, 'sums') == 1 ...   
            || strcmp(NAME_tmp, 'temporalysis') == 1 ...
            || strcmp(NAME_tmp, 'global_count_sum') == 1 ...
            || strcmp(NAME_tmp, 'global_count_sumsquared') == 1 ...
            || strcmp(NAME_tmp, 'radiusCounts') == 1 ...
            || strcmp(NAME_tmp, 'radiusCountsAvg') == 1) ...
            
            continue;
        end
        
        % It was a valid file, so update the NAME and EXT var
        NAME = NAME_tmp;
        
        simnameidx = strfind(NAME, '-');
        simnameidx = str2double(NAME(1:simnameidx-1));

%
%

        createTask(j, @truthradius_onefile, 0, {paramsObj});

        % Submit the job to the scheduler in batches
        if (mod(simCount, TASKS4JOB) == 0)
            submit(j);
            jobCell{jobCount} = j;
            
            if (simCount ~= nCombinations)
                j = createJob(sched); 
                jobCount = jobCount + 1;
            end

        end
        
        
        
%        load(fileName);
%        validFileIdx = validFileIdx + 1;
        
        
        
%         pos = dump.agents;
%         truth = dump.truth;
%         run = dump.run;
%         
%         agentBelong2Radius = zeros(1, nAgents);
%         agentInRadiusFor = zeros(1, nAgents);
%         consensusOnTruth = zeros(1, nIter);
%         consensusHoldsFor = 0;
%         
%         for i = 1:nIter
%         
%             radiusCounts = zeros(1, nRadiusesPlusOne);
%             
%             for a = 1 : nAgents
%                 
%                 found = 0;
%                 
%                 for r = 1 : nRadiuses
%                     radius = RADIUSs(r);
%                     
%                     if (norm(pos(:,a) - truth) <= radius)
%                         found = 1;
%                         break;
%                     end
%                 end
%                 
%                 if (~found)
%                     radius = NOT_IN_RADIUS;
%                     % Incremente r, needs to point to NOT_IN_RADIUS
%                     % in radiusCounts
%                     r = r + 1;
%                 end
%                 
%                 if (agentBelong2Radius(a) ~= radius)
%                     agentBelong2Radius(a) = radius;
%                     agentInRadiusFor(a) = 1;
%                 else
%                     agentInRadiusFor(a) = agentInRadiusFor(a) + 1;                         
%                 end
%                     
%                 if (agentInRadiusFor(a) >= STAY_FOR)
%                     radiusCounts(r) = radiusCounts(r) + 1;
%                 end
%                 
%             end
% 
%             % Is there a consenus on Truth? How long for?
%             if (agentBelong2Radius(TRUTH_RADIUS) > AGENTS_4_CONSENSUS)
%                 consensusHoldsFor = consensusHoldsFor + 1;                
%             else
%                 consensusHoldsFor = 0;
%             end
%             
%             if (consensusHoldsFor > STAY_FOR)
%                 consensusOnTruth(i) = 1;
%             end
%             
%             if (CSV_DUMP)                        
%                 % SAVING ONLY EVERY X ITERATIONS        
%                 if (mod(i, DUMP_RATE) == 0)                                
%                     row_string_final = sprintf(row_string_sprintf, ...
%                         simName, simnameidx, run, i, radiusCount, ...
%                         consensusOnTruth(i));                    
%                     fprintf(fidTruthRadius, '%s\n', row_string_final);   
%                 end
%              
%             end
%             
%             
%             % SUMMING UP AVG statistics
%             globalRadiusCounts(i,:) = globalRadiusCounts(i,:) + radiusCounts;
%             globalRadiusCounts_squared(i,:) = globalRadiusCounts_squared(i,:) + radiusCounts.^2;
%             globalConsensusOnTruth(i) = globalConsensusOnTruth(i) + consensusOnTruth(i)^2;
%             globalConsensusOnTruth_squared(i) = globalConsensusOnTruth_squared(i) + consensusOnTruth(i)^2;
%         end
%         
%         
%         if (PLOTS)
%             
%             % NO PLOTS NOW
%             waitforbuttonpress;
%             
%         end
        
        
    end % File loop
    
    N = validFileIdx; % last value
    df = N - 1;
    
    save([ dumpDir 'radiusCounts' ], ...
                             'globalRadiusCounts', ...
                             'globalRadiusCounts_squared', ...
                             'RADIUSs', ...
                             'STAY_FOR', ...
                             'CONSENSUS_ON_TRUTH_FOR', ...
                             'CONSENSUS_THRESHOLD', ...
                             'AGENTS_4_CONSENSUS', ...
                             'N' ... % number of files
    );
    
    
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
    
    
    
end

%%




