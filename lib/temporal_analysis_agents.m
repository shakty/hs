function temporal_analysis_agents( DUMPDIR, simName, PRECISION, CLU_CUTOFF, CSV_DUMP, DUMP_RATE, PLOTS)

    % Saves the params files, 
    % and other statistiscs on the agents level, clusters are not computed
    
    %% Param
             
  headers_agents = {
        'simname', ...
        'simcount', ...
        'run', ...
        't', ...
        'coverage', ...
        'coverage.cum', ...
        'speed.avg', ...
        'speed.sd', ...
        'move.avg', ...
        'move.sd', ...
        'fromtruth.avg', ...
        'fromtruth.sd' ...
        };
    
    
    headers_agents_avg = {
        'simname', ...
        'simcount', ...
        't', ...
        'coverage.avg', ...
        'coverage.sd', ...
        'coverage.se', ...
        'coverage.ci', ...
        'coverage.cum.avg', ...
        'coverage.cum.sd', ...
        'coverage.cum.se', ...
        'coverage.cum.ci', ...
        'speed.avg', ...
        'speed.sd', ...
        'speed.se', ...
        'speed.ci', ...
        'move.avg', ...
        'move.sd', ...
        'move.se', ...
        'move.ci', ...
        'fromtruth.avg', ...
        'fromtruth.sd', ...
        'fromtruth.se', ...
        'fromtruth.ci' ...
        };
    
   
    headers_params = {
        'simname', ...
        'simcount', ...
        'run', ...
        'timestamp', ...
        't.end', ...
        'dt', ...
        'nagents', ...
        'spacesize',...
        'spacedim', ...
        'alpha', ... % Own velocity
        'R', ...
        'k', ... % Exponent forces
        'A', ... % Attractive force
        'd0', ... 
        'B', ... % Repulsive force
        'd1', ... 
        'tau', ... % Truth strength
        'sigma', ... % Std noise
        'init.vscaling', ...
        'init.nclusters', ...
        'init.clusterradio', ...
        'truth.x', ...
        'truth.y', ...
        'noisetype', ...
        'attrtype', ...
        'attr_on_v', ...
        'seed', ...
        };
    
    
    % CI
    CI_INT = 0.95/2 + 0.5;

    
    % TOTAL CELLS
    TOTAL_CELLS = PRECISION^2;

    
    colnorm = @(X,P) sum(abs(X).^P,1).^(1/P);

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
    % idxsIters = find(mod(1:nIter, DUMP_RATE) == 0);
    % number of files
    nFiles = length(fileIndex);
    
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
    
    % Date and Time
    mytimestamp = datestr ( datevec ( now ), 0 );
    
    if (CSV_DUMP)

        % params (for both)
        paramFileName = [dumpDir 'params.csv'];
        write_csv_headers(paramFileName, headers_params);
        fidParam = fopen(paramFileName,'a');
 
        dataFileName = [dumpDir 'agents_macro.csv'];
        write_csv_headers(dataFileName, headers_agents);
        fidClustersMacro = fopen(dataFileName,'a');
                
        dataFileName = [dumpDir 'agents_macro_avg.csv'];
        write_csv_headers(dataFileName, headers_agents_avg);
        fidClustersMacroAvg = fopen(dataFileName,'a');
       
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
        
        
        load(fileName);
        validFileIdx = validFileIdx + 1;
        
        %dump.parameters
        
        if (CSV_DUMP)
            param_string = csv_format_row_params(simName, simnameidx, dump.run, mytimestamp, dump.parameters, dump.truth);
            % append the param string to the file
            fprintf(fidParam,'%s\n', param_string);
        end
    
        v = dump.agentsv;
        pos = dump.agents;
        
        for i = 1:nIter
            
            % Speed.
            agents_speed = colnorm(v(:,:,i),2);
            mean_agents_speed = mean(agents_speed);
            sd_agents_speed = sd(agents_speed);
            
            % Movement.
            if (i == 1)
                mean_agents_movs = 0;
                sd_agents_movs = 0;
            else
                agents_movs = colnorms(abs(pos(:,:,i) - pos(:,:,i-1)));                
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
            norm_from_truth = norm(pos(:,:,i) - truth);
            mean_agents_fromtruth = mean(norm_from_truth);
            sd_agents_fromtruth = sd(norm_from_truth);
            
            if (CSV_DUMP)
            
            
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
                        'sd_agents_fromtruth', sd_agents_fromtruth ...
                    );
                 
                % 1 Line
                clu_macro_string = csv_format_row_agents_macro(stepData, simName, i);
                fprintf(fidClustersMacro, '%s\n', clu_macro_string);
   
                end
             
            end
            
            
            % SUMMING UP AVG statistics
            global_coverage_sum(i) = global_coverage_sum(i) + avgcoverage;
            global_coverage_sumsquared(i) = global_coverage_sumsquared(i) + avgcoverage^2;
            global_coverage_cum_sum(i) = global_coverage_cum_sum(i) + cumcoverage;
            global_coverage_cum_sumsquared(i) = global_coverage_cum_sumsquared(i) + cumcoverage^2;
            global_speed_sum(i) = global_speed_sum(i) + mean_agents_speed;
            global_speed_sumsquared(i) = global_speed_sumsquared(i) + mean_agents_speed^2;
            global_movs_sum(i) = global_move_sum(i) + mean_agents_movs;
            global_movs_sumsquared(i) = global_move_sumsquared(i) + mean_agents_move^2;
            global_fromtruth_sum(i) = global_fromtruth_sum(i) +  mean_agents_fromtruth;
            global_fromtruth_sumsquared(i) = global_fromtruth_sumsquared(i) + mean_agents_fromtruth^2;
            
    
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
        
        
    end % File loop
    
    N = validFileIdx; % last value
    df = N - 1;
    
    save([ dumpDir 'sums' ], 'global_coverage_sum', ...
                             'global_coverage_sumsquared', ...
                             'global_coverage_cum_sum', ...
                             'global_coverage_cum_sumsquared', ...
                             'global_speed_sum', ...
                             'global_speed_sumsquared', ...
                             'global_movs_sum', ...
                             'global_movs_sumsquared', ...
                             'global_fromtruth_sum', ...
                             'global_fromtruth_sumsquared', ...
                             'N' ... % number of files
    );
    
    
    % Computing global stats  
    
    t_cover_avg = global_coverage_sum / N; 
    t_cover_sd = sqrt(((global_coverage_sumsquared - ((global_coverage_sum).^2 / N))) / df);
    t_cover_se = t_cover_sd / sqrt(N);  
    t_cover_ci = t_cover_se * tquant(CI_INT, df);
    
    t_cover_cum_avg = global_coverage_cum_sum / N; 
    t_cover_cum_sd = sqrt(((global_coverage_cum_sumsquared - ((global_coverage_cum_sum).^2 / N))) / df);
    t_cover_cum_se = t_cover_cum_sd / sqrt(N);  
    t_cover_cum_ci = t_cover_cum_se * tquant(CI_INT, df);
    
    t_speed_avg = global_speed_sum / N; 
    t_speed_sd = sqrt(((global_speed_sumsquared - ((global_speed_sum).^2 / N))) / df);
    t_speed_se = t_speed_sd / sqrt(N);  
    t_speed_ci = t_speed_se * tquant(CI_INT, df);
    
    t_movs_avg = global_movs_sum / N; 
    t_movs_sd = sqrt(((global_movs_sumsquared - ((global_movs_sum).^2 / N))) / df);
    t_movs_se = t_movs_sd / sqrt(N);  
    t_movs_ci = t_movs_se * tquant(CI_INT, df);

    t_fromtruth_avg = global_fromtruth_sum / N; 
    t_fromtruth_sd = sqrt(((global_fromtruth_sumsquared - ((global_fromtruth_sum).^2 / N))) / df);
    t_fromtruth_se = t_fromtruth_sd / sqrt(N);  
    t_fromtruth_ci = t_fromtruth_se * tquant(CI_INT, df);
 
    if (CSV_DUMP)
        
        fclose(fidParam);
        fclose(fidClustersMacro);
        
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
            
            fprintf(fidClustersMacroAvg,'%s\n', clu_macro_avg_string);   
         end

        fclose(fidClustersMacroAvg);
    end
    
    
    
end

%%




