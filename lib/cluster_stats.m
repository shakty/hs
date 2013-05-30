function [ distances, groupCounts, avgGroupDist, avgGroupSpeed, avgGroupMove ] = cluster_stats( T, truth, posagents , v, movs)

    nGroups = max(T);

    avgGroupDist = zeros(nGroups,1);
    groupCounts = zeros(nGroups,1);
    avgGroupSpeed = zeros(nGroups,1);
    avgGroupMove = zeros(nGroups,1);

    distances = zeros(size(posagents,2),1);
    
    for j=1:size(posagents,2)
       distances(j) = norm(posagents(:,j)-truth);
       groupId = T(j);
       groupCounts(groupId) = groupCounts(groupId) + 1;
       avgGroupDist(groupId) = avgGroupDist(groupId) + distances(j);
       avgGroupSpeed(groupId) = avgGroupSpeed(groupId) + v(j);
       avgGroupMove(groupId) = avgGroupMove(groupId) + movs(j);
    end

    avgGroupDist = avgGroupDist ./ groupCounts;
    avgGroupSpeed = avgGroupSpeed ./ groupCounts;
    avgGroupMove = avgGroupMove ./ groupCounts;

end

