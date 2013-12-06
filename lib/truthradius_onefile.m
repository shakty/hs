function [ output_args ] = truthradius_onefile(fileName, simnameidx, RADIUSs, ...
    STAY_FOR, CONSENSUS_ON_TRUTH_FOR, CONSENSUS_THRESHOLD )

    simnameidx = strfind(NAME, '-');
    simnameidx = str2double(NAME(1:simnameidx-1));

    load(fileName);

    pos = dump.agents;
    truth = dump.truth;
    run = dump.run;

    agentBelong2Radius = zeros(1, nAgents);
    agentInRadiusFor = zeros(1, nAgents);
    consensusOnTruth = zeros(1, nIter);
    consensusHoldsFor = 0;

    for i = 1:nIter

        radiusCounts = zeros(1, nRadiusesPlusOne);

        for a = 1 : nAgents

            found = 0;

            for r = 1 : nRadiuses
                radius = RADIUSs(r);

                if (norm(pos(:,a) - truth) <= radius)
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
        if (agentBelong2Radius(TRUTH_RADIUS) > AGENTS_4_CONSENSUS)
            consensusHoldsFor = consensusHoldsFor + 1;                
        else
            consensusHoldsFor = 0;
        end

        if (consensusHoldsFor > STAY_FOR)
            consensusOnTruth(i) = 1;
        end

        if (CSV_DUMP)                        
            % SAVING ONLY EVERY X ITERATIONS        
            if (mod(i, DUMP_RATE) == 0)                                
                row_string_final = sprintf(row_string_sprintf, ...
                    simName, simnameidx, run, i, radiusCount, ...
                    consensusOnTruth(i));                    
                fprintf(fidTruthRadius, '%s\n', row_string_final);   
            end

        end


        % SUMMING UP AVG statistics
        globalRadiusCounts(i,:) = globalRadiusCounts(i,:) + radiusCounts;
        globalRadiusCounts_squared(i,:) = globalRadiusCounts_squared(i,:) + radiusCounts.^2;
        globalConsensusOnTruth(i) = globalConsensusOnTruth(i) + consensusOnTruth(i)^2;
        globalConsensusOnTruth_squared(i) = globalConsensusOnTruth_squared(i) + consensusOnTruth(i)^2;
    end
        
        
        
        

end

