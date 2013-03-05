function [ distances, groupCounts, avgGroupDist ] = each_cluster_from_truth( T, truth, posagents )

    nGroups = max(T);

    avgGroupDist = zeros(nGroups,1);
    groupCounts = zeros(nGroups,1);

    distances = zeros(size(posagents,2),1);

    for j=1:size(posagents,2)
       distances(j) = norm(posagents(:,j)-truth); 
       groupId = T(j);
       groupCounts(groupId) = groupCounts(groupId) + 1;
       avgGroupDist(groupId) = avgGroupDist(groupId) + distances(j);
    end

    avgGroupDist = avgGroupDist ./ groupCounts;
    

end

