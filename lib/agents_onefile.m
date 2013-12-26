function agents_onefile(params)

    % Save statistiscs at the agent level, clusters are not computed.
    
    folderName = params.folderName;
    simName = params.simName;
    fileName = params.fileName;
    PRECISION = params.PRECISION;
    DUMP = params.DUMP;
    DUMP_RATE = params.DUMP_RATE;
    PLOTS = params.PLOTS;
    outDir = params.outDirAgents;
    
    path = [folderName simName fileName];
    load(path);
    
    % Date and Time
    mytimestamp = datestr ( datevec ( now ), 0 );
    
    v = dump.agentsv;
    pos = dump.agents;
    truth = dump.truth;
    run = dump.run;
    simnameidx = dump.sim;
    
    nIter = size(pos,3);
    nAgents = size(pos, 2);
    truth4all = repmat(truth, 1, nAgents);
       
    % TOTAL CELLS
    TOTAL_CELLS = PRECISION^2;
    
    colnorm = @(X,P) sum(abs(X).^P,1).^(1/P);
  
    % GLOBAL variables
  
    global_coverage_sum = zeros(nIter,1);
    global_coverage_sumsquared = zeros(nIter,1);
    global_coverage_cum_sum = zeros(nIter,1);
    global_coverage_cum_sumsquared = zeros(nIter,1);
    global_speed_sum = zeros(nIter,1);
    global_speed_sumsquared = zeros(nIter,1);
    global_movs_sum = zeros(nIter,1);
    global_movs_sumsquared = zeros(nIter,1); 
    global_fromtruth_sum = zeros(nIter,1);
    global_fromtruth_sumsquared = zeros(nIter,1);
    global_pdist_sum =  zeros(nIter,1);
    global_pdist_sumsquared =  zeros(nIter,1);
    
    
    if (DUMP)
        paramFileName = [outDir 'params_' fileName '.csv'];
        fidParam = fopen(paramFileName,'w');
        param_string = csv_format_row_params(simName, simnameidx, run, ...
            mytimestamp, dump.parameters, truth);
        
        % Append the param string to the file
        fprintf(fidParam, '%s\n', param_string);
        fclose(fidParam);
        
        dataFileName = [outDir 'agents_' fileName '.csv'];
        fidAgents = fopen(dataFileName,'w');  
    end

    for i = 1:nIter

        % Speed.
        agents_speed = colnorm(v(:,:,i), 2);
        mean_agents_speed = mean(agents_speed);
        sd_agents_speed = std(agents_speed);

        % Movement.
        if (i == 1)
            mean_agents_movs = 0;
            sd_agents_movs = 0;
        else
            agents_movs = colnorm(pos(:,:,i) - pos(:,:,i-1), 2);              
            mean_agents_movs = mean(agents_movs);
            sd_agents_movs = std(agents_movs);
        end


        % avg share of space
        coverage_matrix = countAgents(pos(:,:,i), PRECISION);
        avgcoverage = nnz(coverage_matrix) / TOTAL_CELLS;

        % cum share of space
        if (i==1)
            cum_coverage_matrix = coverage_matrix;
        else
            cum_coverage_matrix= coverage_matrix + cum_coverage_matrix;
        end
        cumcoverage = nnz(cum_coverage_matrix) / TOTAL_CELLS;

        % avg from truth
        norm_from_truth = colnorm(pos(:,:,i) - truth4all, 2);
        mean_agents_fromtruth = mean(norm_from_truth);
        sd_agents_fromtruth = std(norm_from_truth);
        
        
        % Pair-wise distance between agents.
        pairwise_dist = pdist(pos(:,:,i)', 'euclidean');
        pairwise_dist_mean = mean(pairwise_dist);
        pairwise_dist_sd = std(pairwise_dist);
  
        if (DUMP)

            % SAVING ONLY EVERY X ITERATIONS        
            if (mod(i, DUMP_RATE) == 0)

                stepData = struct(...
                    'simnameidx', simnameidx, ...
                    'run', dump.run, ...
                    'mean_agents_speed', mean_agents_speed, ...
                    'sd_agents_speed', sd_agents_speed, ...
                    'mean_agents_movs', mean_agents_movs, ...
                    'sd_agents_movs', sd_agents_movs, ...
                    'avgcoverage', avgcoverage, ...
                    'cumcoverage', cumcoverage, ...
                    'mean_agents_fromtruth', mean_agents_fromtruth, ... 
                    'sd_agents_fromtruth', sd_agents_fromtruth, ...
                    'pairwise_dist_mean', pairwise_dist_mean, ...
                    'pairwise_dist_sd', pairwise_dist_sd ...
                );

                % 1 Line
                agents_string = csv_format_row_agents_macro(stepData, simName, i);
                fprintf(fidAgents, '%s\n', agents_string);

            end

        end


        % SUMMING UP AVG statistics
        global_coverage_sum(i) = avgcoverage;
        global_coverage_sumsquared(i) = avgcoverage^2;
        global_coverage_cum_sum(i) = cumcoverage;
        global_coverage_cum_sumsquared(i) = cumcoverage^2;
        global_speed_sum(i) = mean_agents_speed;
        global_speed_sumsquared(i) = mean_agents_speed^2;
        global_movs_sum(i) = mean_agents_movs;
        global_movs_sumsquared(i) = mean_agents_movs^2;
        global_fromtruth_sum(i) = mean_agents_fromtruth;
        global_fromtruth_sumsquared(i) = mean_agents_fromtruth^2;
        global_pdist_sum(i) = pairwise_dist_mean;
        global_pdist_sumsquared(i) = pairwise_dist_mean^2;
        

    end


    if (PLOTS)

        close all;
        figure;

        subplot(2,3,1);
        hold on
        plot(1:nIter, avgcoverage, 'r')
        plot(1:nIter, cumcoverage, 'b')
        title('Avg and Cum coverage');
        hold off

        subplot(2,3,2);
        hold on
        plot(1:nIter, mean_agents_speed, 'r')
        plot(1:nIter, sd_agents_speed, 'b')
        title('Mean and Std agents speed');
        hold off

        subplot(2,3,3);
        hold on
        plot(1:nIter, mean_agents_move, 'r')
        plot(1:nIter, sd_agents_move, 'b')
        title('Mean and Std agents move');
        hold off

        subplot(2,3,5);
        hold on
        plot(1:nIter, mean_agents_fromtruth, 'r')
        plot(1:nIter, sd_agents_fromtruth, 'b')
        title('Mean and Std agents from truth');
        hold off


        % Creating a string with the description of the parameters

        paramString = format_sim_params_for_plot_display(simName, NAME, dump.parameters);

        annotation('textbox', [0.7, 0.45, 0, 0], 'string', paramString, ...
            'BackgroundColor', 'white', ...
            'EdgeColor', 'black', ...
            'LineStyle', '-' ...
        );

        waitforbuttonpress;

    end
        
        
    
    if (DUMP)
    
        save([ outDir 'sums_' fileName ], ...
                                    'global_coverage_sum', ...
                                    'global_coverage_sumsquared', ...
                                    'global_coverage_cum_sum', ...
                                    'global_coverage_cum_sumsquared', ...
                                    'global_speed_sum', ...
                                    'global_speed_sumsquared', ...
                                    'global_movs_sum', ...
                                    'global_movs_sumsquared', ...
                                    'global_fromtruth_sum', ...
                                    'global_fromtruth_sumsquared', ...
                                    'global_pdist_sum', ...
                                    'global_pdist_sumsquared' ...
        );
    
    end
   
    
    
    
end

%%




