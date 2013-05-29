%% Saves simulations into properly formatted CSV files

close all;
clear;
clc;

%% Add other directories to path
path(path,'util/'); % Help functions
path(path,'lib/'); % Help functions


CSV_CLU = 0;
CSV_POS = 0;

PLOT_POS = 0;
PLOT_CLU= 0;

colnorm = @(X,P) sum(abs(X).^P,1).^(1/P);

DUMPDIR = 'dump/';
    
simName = 'refactor/refactor-2013-5-29-14-52/';%attr0_av1_nv_seqrnd_Rleft';

dumpDir = [DUMPDIR simName '/'];

simCount = '1-1.mat';

PRECISION = 100;
TOTAL_CELLS = PRECISION^2;

load([dumpDir simCount]);

v = dump.agentsv;
pos = dump.agents;

nIter = size(v,3);
% average velocity of agents at time t
avgspeeds = zeros(1,nIter);
% average distance in space from position at time t-1
avgmovs = zeros(1,nIter);
% average share of space occupied by agents at time t
avgcoverage = zeros(1,nIter);
% cumulative share of space explored by agents at time t
cumcoverage = zeros(1,nIter);
% avg occupation of each cell of the grid at time t
avg_coverage_matrix = zeros(PRECISION, PRECISION, 3);
% cumulative occupation of each cell of the grid at time t
cum_coverage_matrix = zeros(PRECISION, PRECISION, 3);
% cumulative occupation of each cell of the grid at time t
coverage_matrix = zeros(PRECISION, PRECISION, 3);

for i = 1:nIter
    %avg speed
    avgspeeds(i) = mean(colnorm(v(:,:,i),2));
    % avg movement
    if (i>1)
        avgmovs(i) = mean(mean(abs(pos(:,:,i)-pos(:,:,i-1))));
    end
    
    %avg share of space
    coverage_matrix(:,:,i) = countAgents(pos(:,:,i), PRECISION);
    avgcoverage(i) = nnz(coverage_matrix(:,:,i)) / TOTAL_CELLS;
    
    %cum share of space
    if (i==1)
        cum_coverage_matrix(:,:,i) = coverage_matrix(:,:,i);
    else
        cum_coverage_matrix(:,:,i) = coverage_matrix(:,:,i) + cum_coverage_matrix(:,:,i-1);
    end
    cumcoverage(i) = nnz(cum_coverage_matrix(:,:,i)) / TOTAL_CELLS;
    
end

hold on
%plot(1:nIter, avgspeeds, 'r')
%plot(1:nIter, avgmovs, 'b')
plot(1:nIter, avgcoverage, 'g')
plot(1:nIter, cumcoverage, 'g')
hold off


ed de

% Date and Time
mytimestamp = datestr ( datevec ( now ), 0 );

% Load all parameters matrices in one
for i = 1:length(fileIndex)

    append = files(fileIndex(i)).name;
    fileName = [dumpDir, append];

    % We load only .mat
    [PATH,NAME,EXT] = fileparts(fileName);
    if (~strcmpi(EXT,'.mat')) 
        continue;
    end

    simnameidx = strfind(NAME, '-');
    simnameidx = NAME(1:simnameidx-1);

    load(fileName);

    roundsIdx = [1:(1 /dump.parameters.dt):length(dump.agents)];
    rounds = dump.agents(:,:,roundsIdx);

    nRounds = size(rounds,3);

    if (CSV_CLU || CSV_POS)
        param_string = csv_format_row_params(simName, simnameidx, dump.run, mytimestamp, dump.parameters, dump.truth, dump.conv);
        % append the param string to the file
        fprintf(fidParam,'%s\n', param_string);
    end


    Ccount = zeros(1,nRounds);              % Counts the number of clusters
    Csize_mean = zeros(1,nRounds);          % Average cluster size
    Csize_sd = zeros(1,nRounds);          	% Standard deviations of the size of clusters

    Cfromtruth_mean = zeros(1, nRounds);    % Average distance from truth
    Cfromtruth_sd = zeros(1, nRounds);      % Standard deviation of distance from truth

    % Vectors variable length
    % Csize = zeros(1, nRounds);            % Vector of the size of clusters
    % Cfromtruth = zeros(1, nRounds);       % 

    for z=1:nRounds

       agentpos = rounds(:,:,z);

       if (PLOT_POS)
            hold on;
            plot(agentpos(1,:),agentpos(2,:),'rx');
            plot(dump.truth(1),dump.truth(2),'go');
            xlim([0,1])
            ylim([0,1])
            hold off;
            pause(0.1)
       end

       [Z, T, C] = clusterize(agentpos);
       Ccount(1,z) = C;

       [d, Gc, AvgGDist] = each_cluster_from_truth(T, dump.truth, agentpos);

       Csize_mean(1,z) = mean(Gc);
       Csize_sd(1,z) = std(Gc);
       Cfromtruth_mean(1,z) = mean(AvgGDist);
       Cfromtruth_sd(1,z) = std(AvgGDist);


       if (CSV_POS)
            for id = 1:size(agentpos,2)         
                clu_string = csv_format_row_clusters(simName, simnameidx, dump.run, z, id, agentpos(:,id));
                fprintf(fidClusters,'%s\n', clu_string);
            end
       end

       if (CSV_CLU)

            clu_string = csv_format_row_clusters(simName, simnameidx, dump.run, z, ...
                C, Csize_mean(1,z), Csize_sd(1,z), Cfromtruth_mean(1,z), ...
                Cfromtruth_sd(1,z));
            fprintf(fidClusters,'%s\n', clu_string);
       end

    end


    if (PLOT_CLU)
        subplot(2,3,1);
        plot([1:nRounds], Ccount);
        subplot(2,3,2);
        plot([1:nRounds], Csize_mean);
        subplot(2,3,3);
        plot([1:nRounds], Csize_sd);
        subplot(2,3,4);
        plot([1:nRounds], Cfromtruth_mean);
        subplot(2,3,5);
        plot([1:nRounds], Cfromtruth_sd);
        pause(0.5)
    end

end

%%




