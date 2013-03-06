function [ Ccount Csize_mean Csize_sd Cfromtruth_mean Cfromtruth_sd ] = compute_clusters( agentpos )

    [Z, T, Ccount] = clusterize(agentpos);

    [d, Gc, AvgGDist] = each_cluster_from_truth(T, dump.truth, agentpos);

    Csize_mean = mean(Gc);
    Csize_sd = std(Gc);
    Cfromtruth_mean = mean(AvgGDist);
    Cfromtruth_sd = std(AvgGDist);

end

