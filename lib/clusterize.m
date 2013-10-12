function [Z, T, maxC, err] = clusterize(agentpos)

    err = 0;
    X = agentpos';
    Y = pdist(X,'euclidean');
    
    Z = linkage(Y,'average');
    %cophenet(Z,Y);
    %[H,T] = dendrogram(Z,'colorthreshold','default');
    %[H,T] = dendrogram(Z,'colorthreshold',0.1);
    %set(H,'LineWidth',2)
    %T = cluster(Z,'cutoff',1.2)
    T = cluster(Z,'cutoff', 0.1, 'criterion', 'distance');
    maxC = max(T);
    
   
end