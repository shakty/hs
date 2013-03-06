function [ row ] = csv_format_row_clusters( name, nameidx, run, t, count, size_mean, size_sd, fromtruth_mean, fromtruth_sd)
    
    row = [ '"' name '", "' nameidx '", "'  num2str(run) '", ' num2str(t) ] ;
    row = [row, ', ' num2str(count), ', ' num2str(size_mean) ... 
        ', ' num2str(size_sd) ', ' num2str(fromtruth_mean) ...
        ', ' num2str(fromtruth_sd) ];
 
end

