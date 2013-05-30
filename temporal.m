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

% vector of movements of agents at time t
movs = zeros(size(pos,2),1);
% average movement from position at time t-1
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

% Counts the number of clusters
cluster_count = zeros(1,nIter);              
% Average cluster size
mean_cluster_size = zeros(1,nIter);          
% Standard deviations of the size of clusters
sd_cluster_size = zeros(1,nIter);          	
% Average distance from truth
mean_from_truth = zeros(1, nIter);    
% Standard deviation of distance from truth
sd_from_truth = zeros(1, nIter);    


% Following cell arrays containing vector of variable length at each iteration

% Size of the clusters at time t
cluster_sizes = cell(nIter,1);
% Speed of clusters at time t
cluster_speeds = cell(nIter,1);
% Spatial displacement of the cluster compared with time t-1
% It is the average displacement of the agents within in
cluster_movs = cell(nIter,1);
% Distance from truth at time t
cluster_fromtruths = cell(nIter,1);



for i = 1:nIter
    %avg speed
    avgspeeds(i) = mean(colnorm(v(:,:,i),2));
    % avg movement
    if (i>1)
        movs = abs(pos(:,:,i)-pos(:,:,i-1));
        avgmovs(i) = mean(mean(movs));
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
    
    % Z the results of HCLUST 
    % T the cluster of each agent
    % C the number of clusters
    [Z, T, C] = clusterize(pos(:,:,i));
    cluster_count(i) = C;

    [d, Gc, AvgGDist, avgGroupSpeed, avgGroupMove] = cluster_stats(T, dump.truth, pos(:,:,i), v(:,:,i), movs);

    mean_cluster_size(i) = mean(Gc);
    sd_cluster_size(i) = std(Gc);
    mean_from_truth(i) = mean(AvgGDist);
    sd_from_truth(i) = std(AvgGDist);
    
    cluster_sizes{i} = Gc;
    cluster_fromtruths{i} = AvgGDist;
    cluster_speeds{i} = avgGroupSpeed;
    cluster_movs{i} = avgGroupMove;
    
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

    nIter = size(rounds,3);

    if (CSV_CLU || CSV_POS)
        param_string = csv_format_row_params(simName, simnameidx, dump.run, mytimestamp, dump.parameters, dump.truth, dump.conv);
        % append the param string to the file
        fprintf(fidParam,'%s\n', param_string);
    end


    cluster_count = zeros(1,nIter);              % Counts the number of clusters
    mean_cluster_size = zeros(1,nIter);          % Average cluster size
    sd_cluster_size = zeros(1,nIter);          	% Standard deviations of the size of clusters

    mean_from_truth = zeros(1, nIter);    % Average distance from truth
    sd_from_truth = zeros(1, nIter);      % Standard deviation of distance from truth

    % Vectors variable length
    % Csize = zeros(1, nIter);            % Vector of the size of clusters
    % Cfromtruth = zeros(1, nIter);       % 

    for z=1:nIter

       pos = rounds(:,:,z);

       if (PLOT_POS)
            hold on;
            plot(pos(1,:),pos(2,:),'rx');
            plot(dump.truth(1),dump.truth(2),'go');
            xlim([0,1])
            ylim([0,1])
            hold off;
            pause(0.1)
       end

       [Z, T, C] = clusterize(pos);
       cluster_count(1,z) = C;

       [d, Gc, AvgGDist] = each_cluster_from_truth(T, dump.truth, pos);

       mean_cluster_size(1,z) = mean(Gc);
       sd_cluster_size(1,z) = std(Gc);
       mean_from_truth(1,z) = mean(AvgGDist);
       sd_from_truth(1,z) = std(AvgGDist);


       if (CSV_POS)
            for id = 1:size(pos,2)         
                clu_string = csv_format_row_clusters(simName, simnameidx, dump.run, z, id, pos(:,id));
                fprintf(fidClusters,'%s\n', clu_string);
            end
       end

       if (CSV_CLU)

            clu_string = csv_format_row_clusters(simName, simnameidx, dump.run, z, ...
                C, mean_cluster_size(1,z), sd_cluster_size(1,z), mean_from_truth(1,z), ...
                sd_from_truth(1,z));
            fprintf(fidClusters,'%s\n', clu_string);
       end

    end


    if (PLOT_CLU)
        subplot(2,3,1);
        plot([1:nIter], cluster_count);
        subplot(2,3,2);
        plot([1:nIter], mean_cluster_size);
        subplot(2,3,3);
        plot([1:nIter], sd_cluster_size);
        subplot(2,3,4);
        plot([1:nIter], mean_from_truth);
        subplot(2,3,5);
        plot([1:nIter], sd_from_truth);
        pause(0.5)
    end

end

%%




