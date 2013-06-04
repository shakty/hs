function [ row ] = csv_format_row_clusters_micro(stepData, simname, t, cluster_id)
    
    row = [ '"' simname '","' stepData.simnameidx '","' num2str(stepData.run) '","' num2str(t) ] ;
    row = [row, '","' num2str(stepData.clusters_size(cluster_id)), '","' num2str(stepData.clusters_speed(cluster_id)) ];
    row = [row, '","' num2str(stepData.clusters_move(cluster_id)), '","' num2str(stepData.clusters_fromtruth(cluster_id)) '"' ];
    
end

