
% Testing the accuracy of the clustering algorithm

% Add other directories to path
path(path,'../util/'); % Help functions
path(path,'../lib/'); % Help functions

N = 100;
n_agents = [10:10:1000];
nof_cluster = 0;
ideas_space_size = 1;
ideas_space_dim = 2;
clusterTightness = 0;

PLOT = 0;

nSizes = size(n_agents,2);

cs_means_by_nA = zeros(nSizes,1);
cs_stds_by_nA = zeros(nSizes,1);

for j = 1:nSizes
    nA = n_agents(j);
    % Number of clusters found.
    cs = zeros(N,1);
    for i = 1:N
        rng shuffle
        pos = initial_pos_clustered(nof_cluster,clusterTightness, nA,...
                                    ideas_space_size,ideas_space_dim);
        if (PLOT)
            plot(pos(1,:),pos(2,:),'rx');
            pause(0.1);
        end
        [Z, T, C] = clusterize(pos);
        cs(i) = C;
    end
    
    cs_means_by_nA(j) = mean(cs);
    cs_stds_by_nA(j) = std(cs);
end


hold on;
plot(n_agents, cs_means_by_nA, 'rx');
plot(n_agents, cs_stds_by_nA, 'bo');
title('N. of clusters found by n. of agents randomly placed in space (avg over 100 rounds)');
legend('avg', 'std', 'Location', 'SouthEast');
hold off;

% With 1000 repetitions, for 100 agents randomly positioned
% I get on average 41.33 clusters with 2.5 std

% With 1000 repetitions, for 10 agents randomly positioned
% I get on average 8.8 clusters with 0.9 std


% With 1000 repetitions, for 200 agents randomly positioned
% I get on average 52.3 clusters with 2.6 std




