function [ row ] = csv_format_row_clusters_macro( obj , simname, t)
    

    %row = sprintf('"%s","%u","%u","%u","%u","%f","%f","%f","%f","%f","%f","%f","%f","%f","%f","%f"', ... 
    row = sprintf('"%s","%s","%u","%u","%u","%f","%f","%f","%f","%f","%f","%f","%f","%f","%f"', ... 
        simname, ...
        obj.simnameidx, ...
        obj.run, ... 
        t, ...
        obj.cluster_count(t), ...
        obj.avgcoverage(t), ...
        obj.cumcoverage(t), ...
        obj.mean_cluster_speed(t), ...
        obj.sd_cluster_speed(t), ...
        obj.mean_cluster_move(t), ...
        obj.sd_cluster_move(t), ...
        obj.mean_cluster_size(t), ...
        obj.sd_cluster_size(t), ...
        obj.mean_from_truth(t), ...
        obj.sd_from_truth(t) ...
        );
    
    % row = [ '"' simname '","' obj.simnameidx '","'  num2str(obj.run) '","' num2str(t) ] ;
    % row = [row '","' num2str(obj.cluster_count(t)) ];
    % row = [row '","' num2str(obj.avgcoverage(t)) '","'  num2str(obj.cumcoverage(t)) ];
    %row = [row '","' num2str(obj.mean_cluster_speed(t)) '","' num2str(obj.sd_cluster_speed(t)) ];
    %row = [row '","' num2str(obj.mean_cluster_move(t)) '","' num2str(obj.sd_cluster_move(t)) ];
    %row = [row '","'  num2str(obj.mean_cluster_size(t)) '","' num2str(obj.sd_cluster_size(t)) ];
    %row = [row '","'  num2str(obj.mean_from_truth(t)) '","' num2str(obj.sd_from_truth(t)) '"'];
    
end

