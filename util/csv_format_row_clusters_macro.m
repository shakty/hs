function [ row ] = csv_format_row_clusters_macro( obj , simname, t)
    
    row = [ '"' simname '","' obj.simnameidx '","'  num2str(obj.run) '","' num2str(t) ] ;
    row = [row, '","' num2str(obj.cluster_count(t)) '", "aa'];
    row = [row, '","' num2str(obj.avgcoverage(t)) '","'  num2str(obj.cumcoverage(t)) ];
    row = [row, '","' num2str(obj.mean_cluster_speed(t)) '","' num2str(obj.sd_cluster_speed(t)) ];
    row = [row, '","' num2str(obj.mean_cluster_move(t)) '","' num2str(obj.sd_cluster_move(t)) ];
    row = [row, '","'  num2str(obj.mean_cluster_size(t)) '","' num2str(obj.sd_cluster_size(t)) ];
    row = [row, '","'  num2str(obj.mean_from_truth(t)) '","' num2str(obj.sd_from_truth(t)) '"'];
    
end

