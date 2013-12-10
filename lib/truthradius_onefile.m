function truthradius_onefile(params)

    % Counts the number of agents that are distant R from the truth
    % R is every value in the array RADIUSs
    
    % To be considered stable the agent must stay in the same within
    % the same radius for STAY_FOR steps

    folderName = params.folderName;
    simName = params.simName;
    fileName = params.fileName;

    RADIUSs = params.RADIUSs;
    STAY_FOR = params.STAY_FOR;
    CONSENSUS_ON_TRUTH_FOR = params.CONSENSUS_ON_TRUTH_FOR;
    CONSENSUS_THRESHOLD = params.CONSENSUS_THRESHOLD;
    
    DUMP = params.DUMP;
    DUMP_RATE = params.DUMP_RATE;
    outDir = params.outDirRadius;

    PLOTS = params.PLOTS;

    path = [folderName simName fileName];
    load(path);
    
    % Creating variables used to save results to file.
    headers_truthradius = {
        'simname', ...
        'simcount', ...
        'run', ...
        't' ...
    };
   
    row_string_sprintf = '"%s",%u,%u,%u';
    
    % nRadiuses
    nRadiuses = length(RADIUSs);
    hStartFrom = 4;
    for i = 1 : nRadiuses
        radiusStr = ['r_' num2str(RADIUSs(i))];
        headers_truthradius{i + hStartFrom} = radiusStr;
        row_string_sprintf = [row_string_sprintf ',%u'];
    end
    
    % Adding not in radius
    i = i + 1;
    radiusStr = 'r_out';
    headers_truthradius{i + hStartFrom} = radiusStr;
    row_string_sprintf = [row_string_sprintf ',%u'];
    
    % Adding flag consensus on Truth.
    i = i + 1;
    radiusStr = 'consensus';
    headers_truthradius{i + hStartFrom} = radiusStr;
    row_string_sprintf = [row_string_sprintf ',%u']; 
    
    % End creating variables used to save results to file.
    
    pos = dump.agents;
    truth = dump.truth;
    run = dump.run;
    simnameidx = dump.sim;
    
    v = dump.agentsv;
    nIter = size(v,3);
    nAgents = size(v,2);
    AGENTS_4_CONSENSUS = floor(nAgents * CONSENSUS_THRESHOLD);
    
    TRUTH_RADIUS_IDX = 1;
    
    NOT_IN_RADIUS = -1;
    nRadiusesPlusOne = nRadiuses + 1;
    
    agentBelong2Radius = zeros(1, nAgents);
    agentInRadiusFor = zeros(1, nAgents);
    consensusOnTruth = zeros(nIter, 1);
    consensusHoldsFor = 0;
    
    globalRadiusCounts = zeros(nIter, nRadiusesPlusOne);
    globalRadiusCounts_squared = zeros(nIter, nRadiusesPlusOne);
    
    if (DUMP)
        dataFileName = [outDir 'truthradius_' fileName '.csv'];
        write_csv_headers(dataFileName, headers_truthradius);
        fidTruthRadius = fopen(dataFileName, 'w');
    end
    

    for i = 1:nIter

        radiusCounts = zeros(1, nRadiusesPlusOne);

        for a = 1 : nAgents

            found = 0;

            for r = 1 : nRadiuses
                radius = RADIUSs(r);

                if (norm(pos(:, a, i) - truth) <= radius)
                    found = 1;
                    break;
                end
            end

            if (~found)
                radius = NOT_IN_RADIUS;
                % Incremente r, needs to point to NOT_IN_RADIUS
                % in radiusCounts
                r = r + 1;
            end

            if (agentBelong2Radius(a) ~= radius)
                agentBelong2Radius(a) = radius;
                agentInRadiusFor(a) = 1;
            else
                agentInRadiusFor(a) = agentInRadiusFor(a) + 1;                         
            end

            if (agentInRadiusFor(a) >= STAY_FOR)
                radiusCounts(r) = radiusCounts(r) + 1;
            end

        end

        % Is there a consenus on Truth? How long for?
        if (radiusCounts(TRUTH_RADIUS_IDX) > AGENTS_4_CONSENSUS)
            consensusHoldsFor = consensusHoldsFor + 1;                
        else
            consensusHoldsFor = 0;
        end

        if (consensusHoldsFor > CONSENSUS_ON_TRUTH_FOR)
            consensusOnTruth(i) = 1;
        end

        if (DUMP)                        
            % SAVING ONLY EVERY X ITERATIONS        
            if (mod(i, DUMP_RATE) == 0)                                
                row_string_final = sprintf(row_string_sprintf, ...
                    simName, simnameidx, run, i, radiusCounts, ...
                    consensusOnTruth(i));                    
                fprintf(fidTruthRadius, '%s\n', row_string_final);   
            end
        end
        
        
        if (PLOTS)
            %bar(radiusCounts());
            %pause(0.1);
            %waitforbuttonpress
            plot(pos(1,:,i), pos(2,:,i),'rx');
            pause(0.1);            
        end
     
        % SUMMING UP AVG statistics
        globalRadiusCounts(i,:) = globalRadiusCounts(i) + radiusCounts;
        globalRadiusCounts_squared(i,:) = globalRadiusCounts_squared(i) + radiusCounts.^2;
    end
        
    % consensusOnTruth is dicotomous.
    globalConsensusOnTruth = consensusOnTruth;
    globalConsensusOnTruth_squared = consensusOnTruth;
    
    if (DUMP)
        save([ outDir 'sums_' fileName ], ...
                                 'globalRadiusCounts', ...
                                 'globalRadiusCounts_squared', ...
                                 'globalConsensusOnTruth', ...
                                 'globalConsensusOnTruth_squared' ...                                 
        );
    
    end
        
        

end

