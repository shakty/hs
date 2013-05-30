function [ row ] = csv_format_row_clusters_macro( obj )
    
    row = [ '"' obj.name '","' obj.simnameidx '","'  num2str(obj.run) '",' num2str(obj.t) ] ;
    row = [row, '","' num2str(obj.cluster_count) ];
    row = [row, '","' num2str(obj.avgcoverage), '","'  num2str(obj.cumcoverage) ];
    row = [row, '","' num2str(obj.mean_cluster_speed), '","' num2str(obj.sd_cluster_speed) ];
    row = [row, '","' num2str(obj.mean_cluster_move), '","' num2str(obj.sd_cluster_move) ];
    row = [row, '","'  num2str(obj.mean_cluster_size) '","' num2str(obj.sd_cluster_size) ];
    row = [row, '","'  num2str(obj.mean_from_truth) '","' num2str(obj.sd_from_truth) ];
    
end

