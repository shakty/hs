cla
clf

% Testing the accuracy of the clustering algorithm

% Add other directories to path
path(path,'../util/'); % Help functions
path(path,'../lib/'); % Help functions

N = 100;
n_agents = [10];
nof_cluster = 5;
ideas_space_size = 1;
ideas_space_dim = 2;
clusterTightness = [0, 0.01, 0.05, 0.1, 0.15, 0.2];
clusterTightness = 0.2;
PLOT = 1;

nSizes = size(n_agents,2);
nTights = size(clusterTightness,2);
cs_means_by_nA = zeros(nSizes,1);
cs_stds_by_nA = zeros(nSizes,1);

for j = 1:nSizes
    nA = n_agents(j);
    % Number of clusters found.
    cs = zeros(N,1);
    
    for k = 1:nTights
    
        nT = clusterTightness(k);
        
        for i = 1:N
            rng shuffle
            %pos = initial_pos_clustered(nof_cluster,nT, nA,...
            %                            ideas_space_size,ideas_space_dim);
            
            % 101 agents spaced every 0.01 it creates 7 clusters
            % Position in line
            xs = 0:0.01:1;
            pos = [xs; 0.5*ones(1,length(xs))]
            
            if (PLOT)
                clf
                text(pos(1,:),pos(2,:), num2str([1:length(pos)]'));
                %pause(1);
            end
            X = pos';
            hold on
            plot(X(:,1),X(:,2),'rx');
            hold off
            
            Y = pdist(X,'euclidean');
            Ys = squareform(Y)
            Z = linkage(Y,'average');
            %    cophenet(Z,Y);
            %    [H,T] = dendrogram(Z,'colorthreshold','default');
            %    [H,T] = dendrogram(Z,'colorthreshold',0.1);
            %    set(H,'LineWidth',2)
            %T = cluster(Z,'cutoff',1.2)
            T = cluster(Z,'cutoff', 0.1, 'criterion', 'distance')
            maxC = max(T);
            
            cs(i) = C;
        end
        
        cs_means_by_nA(k) = mean(cs);
        cs_stds_by_nA(k) = std(cs);
    end
    
   
end

hold on;
plot(clusterTightness, cs_means_by_nA - 5, 'rx');
plot(clusterTightness, cs_stds_by_nA, 'bo');
title('Difference between clusters detected and placed, by cluster tightness. Agents placed in 5 clusters, with a spread radius of X (avg over 100 rounds)');
legend('avg', 'std', 'Location', 'SouthEast');
hold off;

% With 1000 repetitions, for 100 agents randomly positioned
% I get on average 41.33 clusters with 2.5 std

% With 1000 repetitions, for 10 agents randomly positioned
% I get on average 8.8 clusters with 0.9 std


% With 1000 repetitions, for 200 agents randomly positioned
% I get on average 52.3 clusters with 2.6 std




